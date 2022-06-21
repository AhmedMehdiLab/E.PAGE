#' Read delimited files into tibble
#'
#' @param path path to file
#' @param delim file delimiter, usually \code{","} or \code{"\\t"}
#' @param header whether file contains header
#'
#' @return tibble: file contents
#' @keywords internal
#'
#' @examples \dontrun{
#' path <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
#' file <- import_delim_path(path, ",", TRUE)
#' }
import_delim_path <- function(path, delim, header) {
  type <- readr::cols(.default = readr::col_character())
  file <- readr::read_delim(path, delim, col_types = type, col_names = header,
                            trim_ws = TRUE)

  if (nrow(file) * ncol(file)) return(file)
  stop("File is empty")
}

#' Process tibble of file contents
#'
#' @param file output of \code{\link{import_delim_path}}
#' @param content vector: first and last column of content; default 2 to end
#' @param info column containing information; default none
#'
#' @return tibble: renamed and cleaned file contents
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' path <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
#' file <- import_delim_path(path, ",", TRUE)
#' data <- import_delim_file(file, c(2, 4), 5)
#' }
import_delim_file <- function(file, content, info) {
  if (is.null(content)) lhs <- rhs <- 0
  else {
    lhs <- content[1]
    rhs <- content[2]
  }

  # validate
  . <- NULL
  end <- ncol(file)
  if (lhs < 1 || lhs > end) lhs <- min(2, end)
  if (rhs < 1 || rhs > end) rhs <- end

  # extract
  info <- if (is.null(info) || info < 1 || info > end) character(1)
  else file %>% dplyr::pull(info) %>% replace(is.na(.), "")

  # process
  file %>%
    dplyr::select(name = 1, data_ = dplyr::all_of(lhs:rhs)) %>%
    tibble::add_column(info = info)
}

#' Process tibble of annotation file contents
#'
#' @param anno_file output of \code{\link{import_delim_path}}
#' @param content vector: first and last column of content; default 2 to end
#' @param info column containing information; default none
#'
#' @return glacier-specific imported annotations
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
#' anno_file <- import_delim_path(anno_path, ",", TRUE)
#' anno <- import_annotations_file(anno_file, c(2, 4), 5)
#' }
import_annotations_file <- function(anno_file, content, info) {
  anno_file %>%
    import_delim_file(content, info) %>%
    dplyr::rename(anno_ = dplyr::starts_with("data_"))
}

#' Import delimited annotation file into glacier-specific format
#'
#' @param path path to file
#' @param delim file delimiter, usually \code{","} or \code{"\\t"}
#' @param header whether file contains header
#' @param content vector: first and last column of content; default 2 to end
#' @param info column containing information; default none
#'
#' @return glacier-specific imported annotations
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
#' anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
import_annotations <- function(path, delim, header,
                               content = NULL, info = NULL) {
  path %>%
    import_delim_path(delim, header) %>%
    import_annotations_file(content, info)
}

#' Process tibble of database file contents
#'
#' @param data_file output of \code{\link{import_delim_path}}
#' @param content vector: first and last column of content; default 2 to end
#' @param info column containing information; default none
#'
#' @return glacier-specific imported database
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' data_path <- system.file("extdata", "ex_data.csv", package = "E.PAGE")
#' data_file <- import_delim_path(data_path, ",", FALSE)
#' data <- import_database_file(data_file, c(2, 4), 0)
#' }
import_database_file <- function(data_file, content, info) {
  none <- as.factor("Not assigned")
  proc <- data_file %>%
    import_delim_file(content, info) %>%
    tibble::add_column(category = none, organism = none)

  # extract
  gs_info <- proc %>% dplyr::select(!dplyr::starts_with("data_"))
  gs_genes <- proc %>%
    dplyr::select(dplyr::starts_with("data_")) %>%
    purrr::pmap(c, use.names = F) %>%
    purrr::map(~.[!is.na(.)]) %>%
    purrr::set_names(gs_info$name)

  list(gs_genes = gs_genes, gs_info = gs_info)
}

#' Import delimited database file into glacier-specific format
#'
#' @param path path to file
#' @param delim file delimiter, usually \code{","} or \code{"\\t"}
#' @param header whether file contains header
#' @param content vector: first and last column of content; default 2 to end
#' @param info column containing information; default none
#'
#' @return glacier-specific imported database
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' data_path <- system.file("extdata", "ex_data.csv", package = "E.PAGE")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
import_database <- function(path, delim, header, content = NULL, info = NULL) {
  path %>%
    import_delim_path(delim, header) %>%
    import_database_file(content, info)
}
