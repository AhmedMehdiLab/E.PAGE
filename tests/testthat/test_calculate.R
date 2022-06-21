library(magrittr)
library(tibble)

anno_1 <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
data_1 <- system.file("extdata", "ex_data.csv", package = "E.PAGE")

test_that("calculations can be performed", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  info <- anno[c("name", "info")]

  anno_proc <- process_annotations(anno, info, "file")
  data_proc <- process_database(data)

  input <- process_input_text("")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("FCN1, FTL, CLU")
  expect_equal(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, NULL)$stats$`Odds Ratio`[4],
               5.7842655)

  input <- process_input_text("a 0.1 b 0.3 c 0.8 d e f g 0.9")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("GENE1 GENE2 GENE3 0.2")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, NULL), NA)
})
