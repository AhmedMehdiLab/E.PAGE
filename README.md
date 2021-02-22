# E.PAGE
Environmental Pathways Affecting Gene Expression

![build-macos](https://github.com/AhmedMehdiLab/E-PAGE/actions/workflows/build-macos.yml/badge.svg)

This is a data analysis package allowing the user with a list of genes of interest to find enriched environmental pathways using a curated database.

To run:

```
library(E.PAGE)
input <- process_input_text("GENE1 GENE2 GENE3")
results <- compute(input)
```

Statistically enriched annotations are stored in a `tidyverse` `tibble`, and can be viewed with:

```
results$stats
```
