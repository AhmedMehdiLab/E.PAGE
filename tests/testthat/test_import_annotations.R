library(tibble)

anno_0 <- system.file("extdata", "blank.csv", package = "E.PAGE")
anno_1 <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
data_1 <- system.file("extdata", "ex_data.csv", package = "E.PAGE")

test_that("blank annotation files cannot be imported", {
  expect_error(suppressWarnings(import_annotations(anno_0, ",", F)),
               "File is empty")
})

test_that("annotation files can be imported", {
  expect_equal(
    import_annotations(anno_1, ",", T, c(2, 4), 5),
    tibble(
      name = c("SET_1", "SET_2", "SET_3", "SET_4"),
      anno_1 = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins"),
      anno_2 = c("Viral infection", "Leptin level decrease",
                 "Transcription regulation", "Carcinogen"),
      anno_3 = c("Cell cycle", "Apoptosis", "Cell cycle", NA),
      info = c(
        "Cigarette smoking and viral infections",
        "Genes up-regulated in pancreatic cancer compared to colon cancer",
        "Genes down-regulated in viral infections versus bacterial infections",
        "Genes affected by carcinogens"
      )
    )
  )
  expect_equal(
    import_annotations(anno_1, ",", T, c(2, 4), 0),
    tibble(
      name = c("SET_1", "SET_2", "SET_3", "SET_4"),
      anno_1 = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins"),
      anno_2 = c("Viral infection", "Leptin level decrease",
                 "Transcription regulation", "Carcinogen"),
      anno_3 = c("Cell cycle", "Apoptosis", "Cell cycle", NA),
      info = character(1)
    )
  )
  expect_equal(
    import_annotations(anno_1, ",", T, c(2, 3), 5),
    tibble(
      name = c("SET_1", "SET_2", "SET_3", "SET_4"),
      anno_1 = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins"),
      anno_2 = c("Viral infection", "Leptin level decrease",
                 "Transcription regulation", "Carcinogen"),
      info = c(
        "Cigarette smoking and viral infections",
        "Genes up-regulated in pancreatic cancer compared to colon cancer",
        "Genes down-regulated in viral infections versus bacterial infections",
        "Genes affected by carcinogens"
      )
    )
  )
})
