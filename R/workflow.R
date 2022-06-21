#' Perform computations on imported and processed data
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param anno optional: \code{\link[E.PATH:annotations]{E.PATH::annotations}}
#'   or output of \code{\link{import_annotations}}
#' @param data optional: \code{\link[E.PATH:database]{E.PATH::database}} or
#'   output of \code{\link{import_database}}
#' @param universe number of genes in universe
#' @param info_from optional: \code{"annotation"} or (default) \code{"database"}
#' @param anno_opts optional: \code{"name"} for gene set names, \code{"syms"}
#'   for gene set symbols, \code{"info"} for gene set descriptions,
#'   \code{"auto"} for automatically generated annotations and/or (default)
#'   \code{"file"} for manual annotations
#' @param categories optional: categories to include; default all
#' @param organisms optional: organisms to include; default all
#' @param org_db optional: GO organism database e.g.
#'   \code{\link[org.Hs.eg.db:org.Hs.eg.db]{org.Hs.eg.db::org.Hs.eg.db}}
#' @param save optional: path to save overlap statistics as \code{.csv}
#'
#' @return
#' \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @export
#'
#' @import E.PATH
#' @importFrom magrittr %>%
#' @examples
#' anno <- E.PATH::annotations
#' data <- E.PATH::database
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input)
#' results_custom <- compute(input, anno, data)
compute <- function(input, anno = E.PATH::annotations, data = E.PATH::database,
                    universe = NULL, info_from = "database", anno_opts = "file",
                    categories = FALSE, organisms = FALSE, org_db = NULL,
                    save = NULL) {
  source <- if (info_from == "annotations") anno else data$gs_info
  info <- source %>% dplyr::select("name", "info")

  anno_proc <- process_annotations(anno, info, anno_opts)
  data_proc <- process_database(data, categories, organisms)

  # calculate statistics
  calc <- calculate(input, anno_proc$annos, anno_proc$gs_annos,
                    data_proc$gs_genes, universe, org_db)

  if (!is.null(save)) readr::write_csv(calc$stats, save)
  return(calc)
}
