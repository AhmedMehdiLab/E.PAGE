# E.PAGE
<!-- badges: start -->
[![R-CMD-check](https://github.com/AhmedMehdiLab/E.PAGE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AhmedMehdiLab/E.PAGE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Environmental Pathways Affecting Gene Expression

This is a data analysis package allowing the user with a list of genes of interest to find enriched environmental pathways using a curated database.

## Installation
To install this package, run:

``` r
# install.packages("remotes")
remotes::install_github("AhmedMehdiLab/E.PAGE")
```

## Usage
### Basic usage
To use the package, run:

``` r
library(E.PAGE)

# extract gene list from text input
genes <- "ftl ApoE CTSZ"
input <- process_input_text(genes, capitalize = TRUE)

# alternatively, extract gene list from Seurat object (see documentation for function parameters)
seurat_path <- system.file("extdata", "ex_seurat.rds", package = "E.PAGE")
seurat_obj <- readRDS(seurat_path)
input <- process_input_seurat(seurat_obj, 0)

# compute enriched annotations
results <- compute(input)
```

Statistically enriched annotations are stored in a `tidyverse` `tibble`, and can be viewed with:

``` r
results$stats
```

Roxygen documentation is available for all functions.

### Gene Ontology analysis
Analysis can also be performed on Gene Ontology terms if a database is provided:

```r
library(org.Hs.eg.db)

genes <- "ftl ApoE CTSZ"
input <- process_input_text(genes, capitalize = TRUE)
results <- compute(input, org_db=org.Hs.eg.db)
```

### Custom annotations
This package also supports importing database and annotation files:

```r
# example annotation and database file locations
anno_path <- system.file("extdata/ex_anno.csv", package="E.PAGE")
data_path <- system.file("extdata/ex_data.csv", package="E.PAGE")

# import .csv files (see documentation for function parameters)
anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
data <- import_database(data_path, ",", FALSE, c(2, 4), 0)

# input genes and compute enrichment
genes <- "GENE1 GENE2 GENE3"
input <- process_input_text(genes)
results <- compute(input, anno, data)
```

### Auto-generate annotations
Gene Ontology terms can be automatically generated from database files:

```r
# import and process database file
data_path <- system.file("extdata/ex_data.csv", package="E.PAGE")
data_raw <- import_database(data_path, ",", FALSE, c(2, 4), 0)
data <- process_database(data_raw)

# generate annotations and save to file (see documentation for function parameters)
synthesize_go_anno(data$gs_genes[0:2], limit_universe = FALSE, save = "anno.csv")
```
