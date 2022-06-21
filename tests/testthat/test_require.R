test_that("existing package does not throw error", {
  expect_error(uses("tidyr", stop, "test"), NA)
})

test_that("non-existing package throws error", {
  mockery::stub(uses, "requireNamespace", FALSE)
  expect_error(uses("Seurat", stop, "test"), "test")
})
