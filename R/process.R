# constants
split_1 <- stringr::str_c(sep = "|", " in comparison between ", " upon ",
                          " in comparison of ", " in ", " during ", " after ")
split_2 <- stringr::str_c(sep = "|", " versus ", " vs ", " before ", " after ",
                          " compared to ")
form_lhs <- stringr::str_glue("(?:{split_1})(.*)(?:{split_2})")
form_rhs <- stringr::str_glue("(?:{split_2})(.*?)\\.?$")
utils::globalVariables(c("."))

#' Process annotations
#'
#' Choose annotations to include in annotation list, from gene set names, gene
#' set symbols, gene set descriptions, annotations automatically extracted from
#' gene set descriptions or manual annotations.
#'
#' @param anno output of \code{\link{import_annotations}}
#' @param info tibble: "name" gene set name "info" gene set descriptions
#' @param options \code{"name"} for gene set names, \code{"syms"} for gene set
#'   symbols, \code{"info"} for gene set descriptions, \code{"auto"} for
#'   automatically generated annotations and/or \code{"file"} for manual
#'   annotations
#' @param gs_filter optional: filter for gene set names
#'
#' @return
#' \code{gs_annos} tibble: gene sets and annotations
#'
#' \code{annos} vector: annotations
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
#' anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
#' info <- anno[c("name", "info")]
#'
#' anno_proc <- process_annotations(anno, info, "file")
#' }
process_annotations <- function(anno, info, options, gs_filter="") {
  name <- NULL
  if (gs_filter != "") {
    anno <- anno %>% dplyr::filter(grepl(gs_filter, name))
  }

  anno_proc <- anno %>% dplyr::select("name")
  info <- anno_proc %>%
    dplyr::left_join(info, by = "name") %>%
    dplyr::pull("info")

  # extract annotations
  if ("name" %in% options) anno_proc$anno_name <- anno_proc$name
  if ("syms" %in% options) anno_proc$anno_syms <-
    stringr::str_match(anno_proc$name, "_(.*)")[, 2]
  if ("info" %in% options) anno_proc$anno_info <- info
  if ("auto" %in% options) anno_proc$anno_auto <-
    ifelse(
      stringr::str_detect(info, "up-regulated"),
      stringr::str_match(info, form_lhs)[, 2],
      stringr::str_match(info, form_rhs)[, 2]
    )
  if ("file" %in% options) anno_proc <- anno_proc %>%
    tibble::add_column(dplyr::select(anno, dplyr::starts_with("anno_")))

  # generate annotation list
  annos <- anno_proc %>%
    dplyr::select(dplyr::starts_with("anno_")) %>%
    unlist(use.names = F) %>%
    unique()

  list(gs_annos = anno_proc, annos = annos[!is.na(annos) & annos != ""])
}

#' Process database
#'
#' Choose database category and organism for use in analysis.
#'
#' @param data output of \code{\link{import_database}}
#' @param categories optional: categories to include; default all
#' @param organisms optional: organisms to include; default all
#' @param gs_filter optional: filter for gene set names
#'
#' @return
#' \code{gs_genes} list: names: gene set names vector: genes
#'
#' \code{gs_info} tibble: gene set information
#'
#' \code{genes} vector: list of genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' data_path <- system.file("extdata", "ex_data.csv", package = "E.PAGE")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#' }
process_database <- function(data, categories = FALSE, organisms = FALSE, gs_filter="") {
  # filter categories and organisms
  gs_info <- data$gs_info

  name <- NULL
  if (gs_filter != "") {
    gs_info <- gs_info %>% dplyr::filter(grepl(gs_filter, name))
  }

  if (is.null(categories) || categories != FALSE)
    gs_info <- gs_info %>% dplyr::filter(.data$category %in% categories)
  if (is.null(organisms) || organisms != FALSE)
    gs_info <- gs_info %>% dplyr::filter(.data$organism %in% organisms)

  # extract gene sets and genes
  gs_genes <- data$gs_genes[gs_info$name]
  if (!is.null(gs_genes))
    gs_genes <- gs_genes %>% purrr::modify(toupper)
  genes <- gs_genes %>% unlist(use.names = F) %>% toupper() %>% unique()

  list(gs_genes = gs_genes, gs_info = gs_info, genes = genes)
}

#' Process text input
#'
#' Removes duplicate genes. If multiple values for the same gene are found, only
#' the first value will be kept.
#'
#' @param text character: input with genes and optionally values
#' @param capitalize optional: boolean: force-capitalize all gene names
#'
#' @return tibble: "gene" gene names "value" gene values
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' input <- process_input_text("FCN1 FTL CLU")
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
process_input_text <- function(text, capitalize = FALSE) {
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

  if (capitalize) genes$gene <- toupper(genes$gene)
  return(genes)
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
#' seu_path <- system.file("extdata", "ex_seurat.rds", package = "E.PAGE")
#' seurat <- readRDS(seu_path)
#'
#' input <- process_input_seurat(seurat, 0)
process_input_seurat <- function(seurat, id_1, id_2 = NULL, group = NULL,
                                 cluster = NULL, max_p = 0.05) {
  uses("Seurat", stop, "'Seurat' is required for this feature")
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

#' Generate annotations from data
#'
#' @param data \code{\link[E.PATH:database]{E.PATH::database}} or output of
#'   \code{\link{import_database}}
#' @param mode "GO", "KEGG" or "MeSH"
#' @param limit_universe limit universe to genes contained in \code{data}
#' @param save optional: path to save annotations as \code{.csv}
#'
#' @return tibble: raw annotations
#' @export
#'
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @examples
#' data_path <- system.file("extdata/ex_data.csv", package="E.PAGE")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' auto_anno(data, "GO")
auto_anno <- function(data, mode, limit_universe=FALSE, save = NULL) {
  # set up
  data_proc <- process_database(data)
  univ <- if (limit_universe) data_proc$genes
  ez_univ <- if (limit_universe) AnnotationDbi::mapIds(org.Hs.eg.db, data_proc$genes, "ENTREZID", "SYMBOL")

  find_anno <- function(genes) {
    ez_gene <- AnnotationDbi::mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")

    if (mode == "GO") {
      return(clusterProfiler::enrichGO(ez_gene, "org.Hs.eg.db", ont="ALL", universe=ez_univ))
    } else if (mode == "KEGG") {
      return(clusterProfiler::enrichKEGG(ez_gene, universe=ez_univ))
    } else if (mode == "MeSH") {
      hub <- AnnotationHub::AnnotationHub()
      # hsa <- AnnotationHub::query(hub, c("MeSHDb", "Homo sapiens"))
      mdb <- MeSHDbi::MeSHDb(hub[["AH100340"]]) # Homo sapiens v3

      return(meshes::enrichMeSH(ez_gene, mdb, universe=ez_univ))
    }
  }

  # find annotations
  anno <- data_proc$gs_genes %>%
    purrr::map(toupper) %>%
    purrr::map(find_anno) %>%
    purrr::modify(~ .@result$Description) %>%
    sapply("length<-", max(lengths(.))) %>%
    t() %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "Annotation")

  # store results
  if (!is.null(save)) readr::write_csv(anno, save)
  return(anno)
}

#' Find an annotation's associated gene sets and genes
#'
#' @param annotation annotation to explore
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
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
#' anno <- E.PATH::annotations
#' data <- E.PATH::database
#'
#' anno_assoc <- explore_annotation("Carcinogen", anno$gs_annos, data$gs_genes)
#' }
explore_annotation <- function(annotation, gs_annos, gs_genes, genes = NULL) {
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
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#'
#' @return \code{stats_pre} tibble: overlap statistics (incomplete)
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno <- E.PATH::annotations
#' data <- E.PATH::database
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' calc_pre <- calculate_pre(input, anno$annos, anno$gs_annos, data$gs_genes)
#' }
calculate_pre <- function(input, annos, gs_annos, gs_genes) {
  stat <- tibble::tibble(name = annos, n_sets = 0L, n_gene = 0L, n_hits = 0L,
                         p_hits = 0, g_hits = "")
  hits <- list()

  # iterate over annotations
  for (i in seq_along(annos)) {
    # get related genes and find overlap
    overlap <- explore_annotation(annos[i], gs_annos, gs_genes)
    matches <- intersect(overlap$genes, input$gene)

    # store information
    stat$n_sets[i] <- length(overlap$names)
    stat$n_gene[i] <- length(overlap$genes)
    stat$n_hits[i] <- length(matches)
    stat$p_hits[i] <- length(matches) / nrow(input) * 100
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
#' anno <- E.PATH::annotations
#' data <- E.PATH::database
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#'
#' calc_pre <- calculate_pre(input, anno$annos, anno$gs_annos, data$gs_genes)
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
    `% match` = .data$p_hits,
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
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#' @param univ_size optional: number of genes in universe; default calculated
#'   from \code{gs_genes}
#' @param org_db optional: GO organism database e.g.
#'   \code{\link[org.Hs.eg.db:org.Hs.eg.db]{org.Hs.eg.db::org.Hs.eg.db}}
#'
#' @return \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#'
#' \code{go} \code{NULL} or output of
#'   \code{\link[clusterProfiler:enrichGO]{clusterProfiler::enrichGO}}
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno <- E.PATH::annotations
#' data <- E.PATH::database
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' calc <- calculate(input, anno$annos, anno$gs_annos, data$gs_genes)
#' }
calculate <- function(input, annos, gs_annos, gs_genes, univ_size = NULL,
                      org_db = NULL) {
  # extract background gene list
  background <- gs_genes %>% unlist(use.names = F)

  # calculate size of universe
  if (is.null(univ_size)) univ_size <- background %>%
      c(input$gene) %>%
      toupper() %>%
      unique() %>%
      length()

  calc_pre <- calculate_pre(input, annos, gs_annos, gs_genes)
  calc <- calculate_post(calc_pre$stats_pre, nrow(input), univ_size)

  go <- if (is.null(org_db)) NULL
  else clusterProfiler::enrichGO(
    toupper(input$gene), org_db, "SYMBOL", "ALL", universe = toupper(background)
  )

  list(stats = calc, matches = calc_pre$matches, go = go)
}
