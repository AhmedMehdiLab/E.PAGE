seu_path <- system.file("extdata", "ex_seurat.rds", package = "E.PAGE")
seurat <- readRDS(seu_path)

test_that("differentially expressed genes can be found", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seurat, 0)$gene,
               c("FTL", "APOE", "CTSZ", "C1QA"))
})

test_that("correctly handles case when two identical clusters specified", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seu, 0, 0)$gene, character())
})
