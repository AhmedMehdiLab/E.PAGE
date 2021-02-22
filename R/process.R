#' Process text input
#'
#' Removes duplicate genes. If multiple values for the same gene are found, only
#' the first value will be kept.
#'
#' @param text character: input with genes and optionally values
#'
#' @return tibble: "gene" gene names "value" gene values
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' input <- process_input_text("FCN1 FTL CLU")
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
process_input_text <- function(text) {
  tokens <- text %>%
    stringr::str_split("[ \t\r\n,;]") %>%
    unlist() %>%
    purrr::discard(~. == "")
  values <- suppressWarnings(as.numeric(tokens))

  # process
  tokens[!is.na(values)] <- NA
  values <- values[-1] %>% c(NA)

  # store results
  genes <- tibble::tibble(gene = tokens, value = values) %>%
    tidyr::drop_na(.data$gene) %>%
    dplyr::distinct(.data$gene, .keep_all = T)
}

#' Extract differentially expressed genes from Seurat object
#'
#' Finds differentially expressed genes, records adjusted P-value and filters
#' for values less than \code{max_p}.
#'
#' @param seurat Seurat object
#' @param id_1 first identity
#' @param id_2 optional: second identity; default all others
#' @param group optional: subgroup of cluster
#' @param cluster optional: cluster selected
#' @param max_p P-value cutoff, only genes with P-values less than this will be
#'   returned
#'
#' @return tibble: "gene" gene names "value" gene values
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' seurat <- E.PAGE:::seurat
#' input <- process_input_seurat(seurat, 0)
process_input_seurat <- function(seurat, id_1, id_2 = NULL, group = NULL,
                                 cluster = NULL, max_p = 0.05) {
  if (!is.null(id_2) && id_1 == id_2)
    return(tibble::tibble(gene = character(), value = numeric()))

  seurat %>%
    Seurat::FindMarkers(ident.1 = id_1, ident.2 = id_2, group.by = group,
                        subset.ident = cluster) %>%
    tibble::rownames_to_column("gene") %>%
    tibble::tibble() %>%
    dplyr::select("gene", value = "p_val_adj") %>%
    dplyr::filter(.data$value <= max_p)
}

#' Find an annotation's associated gene sets and genes
#'
#' @param annotation annotation to explore
#' @param genes optional: genes to match, or (default) all
#'
#' @return
#' \code{"names"} vector: gene set names
#'
#' \code{"genes"} vector: genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_assoc <- explore_annotation("Carcinogen")
#' }
explore_annotation <- function(annotation, genes = NULL) {
  gs_annos <- anno$gs_annos
  gs_genes <- data$gs_genes

  index <- (gs_annos == annotation) %>% rowSums(na.rm = T) %>% as.logical()
  match <- gs_genes[gs_annos$name[index]]
  if (!is.null(genes))
    match <- match %>% purrr::map(~intersect(., genes)) %>% purrr::compact()

  names <- names(match)
  genes <- match %>% unlist(use.names = F) %>% unique()
  list(names = names, genes = if (is.null(genes)) character() else genes)
}

#' Begin calculating overlap statistics
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#'
#' @return \code{stats_pre} tibble: overlap statistics (incomplete)
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' calc_pre <- calculate_pre(input)
#' }
calculate_pre <- function(input) {
  annos <- anno$annos
  gs_annos <- anno$gs_annos
  gs_genes <- data$gs_genes

  stat <- tibble::tibble(name = annos, n_sets = 0L, n_gene = 0L, n_hits = 0L,
                         g_hits = "")
  hits <- list()

  # iterate over annotations
  for (i in seq_along(annos)) {
    # get related genes and find overlap
    overlap <- explore_annotation(annos[i])
    matches <- intersect(overlap$genes, input$gene)

    # store information
    stat$n_sets[i] <- length(overlap$names)
    stat$n_gene[i] <- length(overlap$genes)
    stat$n_hits[i] <- length(matches)
    stat$g_hits[i] <- matches %>% stringr::str_c(collapse = ", ")
    hits[[annos[i]]] <- matches
  }

  list(stats_pre = stat, matches = hits)
}

#' Finish calculating overlap statistics
#'
#' @param stats_pre value from \code{\link{calculate_pre}}
#' @param input_size number of genes in input
#' @param universe number of genes in universe
#'
#' @return tibble: overlap statistics
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' calc_pre <- calculate_pre(input)
#' calc <- calculate_post(calc_pre$stats_pre, nrow(input), 10000)
#' }
calculate_post <- function(stats_pre, input_size, universe) {
  stat <- stats_pre %>% tibble::add_column(pvalue = 0, odds_r = 0)

  # calculate statistics
  for (i in seq_len(nrow(stat))) {
    # fisher's exact test: [1, ] belongs to annotation [, 1] entered in list
    data <- matrix(nrow = 2, ncol = 2)
    data[1, 1] <- stat$n_hits[i]
    data[1, 2] <- stat$n_gene[i] - data[1, 1]
    data[2, 1] <- input_size - data[1, 1]
    data[2, 2] <- universe - data[1, 1] - data[1, 2] - data[2, 1]

    # assign statistics
    test <- data %>% stats::fisher.test(alternative = "greater")
    stat$pvalue[i] <- test$p.value
    stat$odds_r[i] <- test$estimate[[1]]
  }

  # post-processing
  stat$enrich <- (stat$n_hits / input_size) / (stat$n_gene / universe)
  stat$adj_pv <- stat$pvalue %>% stats::p.adjust(method = "fdr")
  stat$adj_fe <- stat$enrich / -log(stat$adj_pv)

  stat %>% dplyr::select(
    Annotation = .data$name,
    `# gene sets` = .data$n_sets,
    `# genes` = .data$n_gene,
    `# matches` = .data$n_hits,
    `P-value` = .data$pvalue,
    `Adjusted P-value` = .data$adj_pv,
    `Odds Ratio` = .data$odds_r,
    `Fold Enrichment` = .data$enrich,
    `Adjusted Fold Enrichment` = .data$adj_fe,
    Matches = .data$g_hits
  )
}

#' Calculate overlap statistics
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param universe optional: number of genes in universe; default calculate from
#'   \code{gs_genes}
#'
#' @return \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' calc <- calculate(input)
#' }
calculate <- function(input, universe = NULL) {
  gs_genes <- data$gs_genes

  if (is.null(universe)) universe <- gs_genes %>%
      unlist(use.names = F) %>%
      c(input$gene) %>%
      unique() %>%
      length()

  calc_pre <- calculate_pre(input)
  calc <- calculate_post(calc_pre$stats_pre, nrow(input), universe)

  list(stats = calc, matches = calc_pre$matches)
}
