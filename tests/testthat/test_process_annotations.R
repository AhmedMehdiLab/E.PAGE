library(dplyr)
library(magrittr)

anno_1 <- system.file("extdata", "ex_anno.csv", package = "E.PAGE")
anno_info_ex <- c(
  "Cigarette smoking and viral infections",
  "Genes up-regulated in pancreatic cancer compared to colon cancer",
  "Genes down-regulated in viral infections versus bacterial infections",
  "Genes affected by carcinogens"
)

test_that("annotations can be processed", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  info <- anno[c("name", "info")]

  expect_equal(
    process_annotations(anno, info, ""),
    list(gs_annos = anno %>% select("name"), annos = NULL)
  )
  expect_equal(
    process_annotations(anno, info, "name"),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(anno_name = c("SET_1", "SET_2", "SET_3", "SET_4")),
      annos = c("SET_1", "SET_2", "SET_3", "SET_4")
    )
  )
  expect_equal(
    process_annotations(anno, info, "syms"),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(anno_syms = c("1", "2", "3", "4")),
      annos = c("1", "2", "3", "4")
    )
  )
  expect_equal(
    process_annotations(anno, info, "info"),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(anno_info = anno_info_ex),
      annos = anno_info_ex
    )
  )
  expect_equal(
    process_annotations(anno, info, "auto"),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(
          anno_auto = c(NA, "pancreatic cancer", "bacterial infections", NA)
        ),
      annos = c("pancreatic cancer", "bacterial infections")
    )
  )
  expect_equal(
    process_annotations(anno, info, "file"),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(
          anno_1 = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins"),
          anno_2 = c("Viral infection", "Leptin level decrease",
                     "Transcription regulation", "Carcinogen"),
          anno_3 = c("Cell cycle", "Apoptosis", "Cell cycle", NA)
        ),
      annos = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins",
                "Viral infection", "Leptin level decrease",
                "Transcription regulation", "Carcinogen", "Cell cycle")
    )
  )
  expect_equal(
    process_annotations(anno, info,
                        c("name", "syms", "info", "auto", "file")),
    list(
      gs_annos = anno %>%
        select("name") %>%
        add_column(
          anno_name = c("SET_1", "SET_2", "SET_3", "SET_4"),
          anno_syms = c("1", "2", "3", "4"),
          anno_info = anno_info_ex,
          anno_auto = c(NA, "pancreatic cancer", "bacterial infections", NA),
          anno_1 = c("Asthma", "Insulin release", "Apoptosis", "Neurotoxins"),
          anno_2 = c("Viral infection", "Leptin level decrease",
                     "Transcription regulation", "Carcinogen"),
          anno_3 = c("Cell cycle", "Apoptosis", "Cell cycle", NA)
        ),
      annos = c(
        "SET_1", "SET_2", "SET_3", "SET_4", "1", "2", "3", "4", anno_info_ex,
        "pancreatic cancer", "bacterial infections", "Asthma",
        "Insulin release", "Apoptosis", "Neurotoxins", "Viral infection",
        "Leptin level decrease", "Transcription regulation", "Carcinogen",
        "Cell cycle"
      )
    )
  )
})
