library(tibble)

data_0 <- system.file("extdata", "blank.csv", package = "E.PAGE")
data_1 <- system.file("extdata", "ex_data.csv", package = "E.PAGE")

test_that("blank databases cannot be imported", {
  expect_error(suppressWarnings(import_database(data_0, ",", F)),
               "File is empty")
})

test_that("databases can be imported", {
  expect_equal(
    import_database(data_1, ",", F, c(2, 4), 0),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL", "CLU"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("Bim", "CLU", "FoxO3a"),
                      SET_4 = c("FTL", "SOX2", "FCN1")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = character(4),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      )
    )
  )
  expect_equal(
    import_database(data_1, ",", F),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL", "CLU"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("Bim", "CLU", "FoxO3a"),
                      SET_4 = c("FTL", "SOX2", "FCN1")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = character(4),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      )
    )
  )
  expect_equal(
    import_database(data_1, ",", F, c(2, 4), 1),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL", "CLU"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("Bim", "CLU", "FoxO3a"),
                      SET_4 = c("FTL", "SOX2", "FCN1")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = c("SET_1", "SET_2", "SET_3", "SET_4"),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      )
    )
  )
  expect_equal(
    import_database(data_1, ",", F, c(2, 3), 0),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("Bim", "CLU"),
                      SET_4 = c("FTL", "SOX2")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = character(1),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      )
    )
  )
})
