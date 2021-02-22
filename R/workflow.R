#' Perform computations on imported and processed data
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param universe number of genes in universe
#' @param save optional: path to save overlap statistics as \code{.csv}
#'
#' @return
#' \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input)
compute <- function(input, universe = NULL, save = NULL) {
  calc <- calculate(input, universe)
  if (!is.null(save)) readr::write_csv(calc$stats, save)
  return(calc)
}
