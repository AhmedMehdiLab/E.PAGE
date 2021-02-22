#' Plot gene-annotation matches as heatmap
#'
#' Heatmap will plot genes on the X-axis against annotations on the Y-axis.
#' Matches will be shown with colored squares, and non-matches will be
#' transparent. If the \code{value} is "Gene Value", "Odds Ratio", "Fold
#' Enrichment" or "Adjusted Fold Enrichment", if the value is NaN or an
#' infinity, these will also be transparent.
#'
#' @param matches value from \code{\link{compute}}
#' @param value \code{"Gene Value"} for values from \code{input} or one of
#'   \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"} for values
#'   from \code{stats}
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param stats value from \code{\link{compute}}
#' @param xpos x axis position, "top" or "bottom"
#'
#' @return ggplot2: heatmap of overlap
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_raster
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input)
#' over <- plot_overlap(results$matches, "Gene Value", input, results$stats)
plot_overlap <- function(matches, value, input, stats, xpos = "bottom") {
  . <- NULL
  genes <- input$gene %>% factor(., levels = .)
  annos <- stats$Annotation %>% factor(., levels = .) %>% forcats::fct_rev()

  # helper function to find gene values
  get_value <- function(g, a) {
    if (value == "Gene Value") input %>%
      dplyr::filter(.data$gene == g) %>%
      dplyr::pull(.data$value)
    else stats %>%
      dplyr::filter(.data$Annotation == a) %>%
      dplyr::pull(!!value)
  }

  # construct data grid
  data <- expand.grid(Gene = genes, Annotation = annos, Value = NA_real_)
  for (i in seq_len(nrow(data))) {
    anno <- data$Annotation[i] %>% as.character()
    gene <- data$Gene[i] %>% as.character()

    if (gene %in% matches[[anno]]) data$Value[i] <- get_value(gene, anno)
  }
  colnames(data)[3] <- value

  # plot data
  ggplot(data, aes(.data$Gene, .data$Annotation, fill = .data[[value]])) +
    geom_raster() + if (nrow(data)) ggplot2::scale_x_discrete(position = xpos)
}

#' Plot overlap statistics as a bar graph
#'
#' @param stats value from \code{\link{compute}}
#' @param value \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}
#' @param color \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}
#' @param sort_y whether to sort annotations by \code{value}
#'
#' @return ggplot2: bar chart of statistics
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @examples
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input)
#' stat <- plot_stats(results$stats, 'Fold Enrichment', 'Adjusted P-value')
plot_stats <- function(stats, value, color, sort_y = FALSE) {
  # prepare axes
  . <- NULL
  value <- rlang::sym(value)
  color <- rlang::sym(color)

  # order Y axis
  if (nrow(stats) != 0 && sort_y) stats %<>% dplyr::arrange(dplyr::desc(!!value))
  stats$Annotation %<>% factor(., levels = .) %>% forcats::fct_rev()

  # plot data
  ggplot(stats, aes(!!value, .data$Annotation, fill = !!color)) + geom_col()
}
