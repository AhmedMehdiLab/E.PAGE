library(org.Hs.eg.db)

test_that("Basic usage works", {
  genes <- "ftl ApoE CTSZ"
  input <- process_input_text(genes, capitalize = TRUE)

  seurat_path <- system.file("extdata", "ex_seurat.rds", package = "E.PAGE")
  seurat_obj <- readRDS(seurat_path)
  input <- process_input_seurat(seurat_obj, 0)

  results <- compute(input)
  expect_true(!is.null(results$stats))
})

test_that("Gene Ontology analysis works", {
  genes <- "ftl ApoE CTSZ"
  input <- process_input_text(genes, capitalize = TRUE)
  results <- compute(input, org_db=org.Hs.eg.db)
  expect_true(!is.null(results$go))
})

test_that("Custom annotations work", {
  anno_path <- system.file("extdata/ex_anno.csv", package="E.PAGE")
  data_path <- system.file("extdata/ex_data.csv", package="E.PAGE")

  anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
  data <- import_database(data_path, ",", FALSE, c(2, 4), 0)

  genes <- "GENE1 GENE2 GENE3"
  input <- process_input_text(genes)
  results <- compute(input, anno, data)
  expect_true(!is.null(results))
})

test_that("Auto-generate annotations works", {
  data_path <- system.file("extdata/ex_data.csv", package="E.PAGE")
  data_raw <- import_database(data_path, ",", FALSE, c(2, 4), 0)
  data <- process_database(data_raw)

  go <- synthesize_go_anno(data$gs_genes[0], limit_universe = FALSE, save = "anno.csv")
  expect_true(!is.null(go))
})
