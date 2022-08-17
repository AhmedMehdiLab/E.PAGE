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
#' @param gs_filter optional: filter for gene set names
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
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input)
compute <- function(input, anno = E.PATH::annotations, data = E.PATH::database,
                    universe = NULL, info_from = "database", anno_opts = "file",
                    categories = FALSE, organisms = FALSE, org_db = NULL,
                    save = NULL, gs_filter="") {
  source <- if (info_from == "annotations") anno else data$gs_info
  info <- source %>% dplyr::select("name", "info")

  anno_proc <- process_annotations(anno, info, anno_opts, gs_filter)
  data_proc <- process_database(data, categories, organisms, gs_filter)

  # calculate statistics
  calc <- calculate(input, anno_proc$annos, anno_proc$gs_annos,
                    data_proc$gs_genes, universe, org_db)

  if (!is.null(save)) readr::write_csv(calc$stats, save)
  return(calc)
}

#' Compute overlap enrichment on E.PATH gene set categories
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
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
#'
#' @return list: categories and output of \code{\link{compute}}
#' @export
#'
#' @examples
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute_sets(input)
compute_sets <- function(input, universe = NULL, info_from = "database",
                         anno_opts = "file", categories = FALSE,
                         organisms = FALSE, org_db = NULL) {
  return(list(
    CSGE = compute(input, universe=universe, info_from=info_from, anno_opts=anno_opts,
                   categories=categories, organisms=organisms, org_db=org_db, gs_filter="CSGE"),
    IGE = compute(input, universe=universe, info_from=info_from, anno_opts=anno_opts,
                   categories=categories, organisms=organisms, org_db=org_db, gs_filter="IGE"),
    DGE = compute(input, universe=universe, info_from=info_from, anno_opts=anno_opts,
                   categories=categories, organisms=organisms, org_db=org_db, gs_filter="DGE"),
    TGE = compute(input, universe=universe, info_from=info_from, anno_opts=anno_opts,
                   categories=categories, organisms=organisms, org_db=org_db, gs_filter="TGE")
  ))
}
