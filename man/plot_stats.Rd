% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R
\name{plot_stats}
\alias{plot_stats}
\title{Plot overlap statistics as a bar graph}
\usage{
plot_stats(stats, value, color, sort_y = FALSE)
}
\arguments{
\item{stats}{value from \code{\link{compute}}}

\item{value}{\code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
\code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
\code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}}

\item{color}{\code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
\code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
\code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}}

\item{sort_y}{whether to sort annotations by \code{value}}
}
\value{
ggplot2: bar chart of statistics
}
\description{
Plot overlap statistics as a bar graph
}
\examples{
anno <- E.PATH::annotations
data <- E.PATH::database

input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
results <- compute(input, anno, data)
stat <- plot_stats(results$stats, 'Fold Enrichment', 'Adjusted P-value')
}
