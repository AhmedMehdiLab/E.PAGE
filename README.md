# E.PAGE
<!-- badges: start -->
[![R-CMD-check](https://github.com/AhmedMehdiLab/E.PAGE/workflows/R-CMD-check/badge.svg)](https://github.com/AhmedMehdiLab/E.PAGE/actions)
<!-- badges: end -->

Environmental Pathways Affecting Gene Expression

This is a data analysis package allowing the user with a list of genes of interest to find enriched environmental pathways using a curated database.

## Installation
To install this package, run:

``` r
remotes::install_github("AhmedMehdiLab/E.PAGE")
```

## Usage
To use the package, run:

``` r
library(E.PAGE)
genes <- "GENE1 GENE2 GENE3"       # a gene list separated by spaces or commas
input <- process_input_text(genes) 
results <- compute(input)
```

Statistically enriched annotations are stored in a `tidyverse` `tibble`, and can be viewed with:

```
results$stats
```

Roxygen documentation is available for all functions.
