% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{process_database}
\alias{process_database}
\title{Process database}
\usage{
process_database(data, categories = FALSE, organisms = FALSE, gs_filter = "")
}
\arguments{
\item{data}{output of \code{\link{import_database}}}

\item{categories}{optional: categories to include; default all}

\item{organisms}{optional: organisms to include; default all}

\item{gs_filter}{optional: filter for gene set names}
}
\value{
\code{gs_genes} list: names: gene set names vector: genes

\code{gs_info} tibble: gene set information

\code{genes} vector: list of genes
}
\description{
Choose database category and organism for use in analysis.
}
\examples{
\dontrun{
data_path <- system.file("extdata", "ex_data.csv", package = "E.PAGE")
data <- import_database(data_path, ",", FALSE, c(2, 4), 0)

data_proc <- process_database(data, "Not assigned", "Not assigned")
}
}
\keyword{internal}
