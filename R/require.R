#' Checks that package can be used, else returns error
#'
#' @param packages character vector: packages required
#' @param error error function
#' @param message error message
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' uses("stringr", stop, "Prerequisite not available")
#' }
uses <- function(packages, error, message) {
  for (package in packages)
    if (!requireNamespace(package, quietly = TRUE))
      error(message)
}