% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{explore_annotation}
\alias{explore_annotation}
\title{Find an annotation's associated gene sets and genes}
\usage{
explore_annotation(annotation, gs_annos, gs_genes, genes = NULL)
}
\arguments{
\item{annotation}{annotation to explore}

\item{gs_annos}{value from \code{\link{process_annotations}}}

\item{gs_genes}{value from \code{\link{process_database}}}

\item{genes}{optional: genes to match, or (default) all}
}
\value{
\code{"names"} vector: gene set names

\code{"genes"} vector: genes
}
\description{
Find an annotation's associated gene sets and genes
}
\examples{
\dontrun{
anno <- E.PATH::annotations
data <- E.PATH::database

anno_assoc <- explore_annotation("Carcinogen", anno$gs_annos, data$gs_genes)
}
}
\keyword{internal}
