library(magrittr)
library(purrr)
library(tibble)

data_1 <- system.file("extdata", "ex_data.csv", package = "E.PAGE")

test_that("complex databases can be processed", {
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  none <- factor("Not assigned")

  expect_equal(
    process_database(data, NULL, NULL),
    list(
      gs_genes = list() %>% set_names(),
      gs_info = tibble(name = character(), info = character(), category = none,
                       organism = none),
      genes = character()
    )
  )
  expect_equal(
    process_database(data, "Not assigned", "Not assigned"),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL", "CLU"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("BIM", "CLU", "FOXO3A"),
                      SET_4 = c("FTL", "SOX2", "FCN1")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = character(4),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      ),
      genes = c("FCN1", "FTL", "CLU", "PDK4", "PFKFB3", "BIM", "FOXO3A", "SOX2")
    )
  )
})
