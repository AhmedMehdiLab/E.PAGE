library(magrittr)
library(tibble)

base <- tibble(gene = character(), value = numeric())

test_that("works for well-formed single gene input", {
  expect_equal(process_input_text(""), base)
  expect_equal(process_input_text("a"), base %>% add_row(gene = "a"))
  expect_equal(process_input_text("0.1"), base)
  expect_equal(process_input_text("a 0.1"),
               base %>% add_row(gene = "a", value = 0.1))

  expect_equal(process_input_text(" "), base)
  expect_equal(process_input_text("a "), base %>% add_row(gene = "a"))
  expect_equal(process_input_text("0.1 "), base)
  expect_equal(process_input_text("a 0.1 "),
               base %>% add_row(gene = "a", value = 0.1))
})

test_that("works for well-formed multiple gene input", {
  expect_equal(process_input_text("a b"),
               base %>% add_row(gene = "a") %>% add_row(gene = "b"))
  expect_equal(process_input_text("a b c"), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b") %>%
                 add_row(gene = "c"))
  expect_equal(process_input_text("a 0.1 b 0.2"), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b", value = 0.2))
  expect_equal(process_input_text("a 0.1 b 0.2 c 0.1"), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b", value = 0.2) %>%
                 add_row(gene = "c", value = 0.1))

  expect_equal(process_input_text("a b "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b"))
  expect_equal(process_input_text("a b c "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b") %>%
                 add_row(gene = "c"))
  expect_equal(process_input_text("a 0.1 b 0.1"), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b", value = 0.1))
  expect_equal(process_input_text("a 0.1 b 0.1 c 0.1 "), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b", value = 0.1) %>%
                 add_row(gene = "c", value = 0.1))
})

test_that("eliminates duplicate genes", {
  expect_equal(process_input_text(" "), base)
  expect_equal(process_input_text("  "), base)
  expect_equal(process_input_text("\n\r\t \n\r\t "), base)

  expect_equal(process_input_text("a a"), base %>% add_row(gene = "a"))
  expect_equal(process_input_text("a b a "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b"))
  expect_equal(process_input_text("a a b b "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b"))
  expect_equal(process_input_text("a b a b a b "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b"))
  expect_equal(process_input_text("a b c b a a a b b c c c a b c "), base %>%
                 add_row(gene = "a") %>%
                 add_row(gene = "b") %>%
                 add_row(gene = "c"))
})

test_that("eliminates duplicate values", {
  expect_equal(process_input_text("a 0.1 0.2 0.3"), base %>%
                 add_row(gene = "a", value = 0.1))
  expect_equal(process_input_text("a 0.1 b 0.2 a 0.2 b 0.4 a 0.2 0.4"), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b", value = 0.2))
  expect_equal(process_input_text("a 0.1 b c 0.1 b 0.1 c 0.4 a 0.8"), base %>%
                 add_row(gene = "a", value = 0.1) %>%
                 add_row(gene = "b") %>%
                 add_row(gene = "c", value = 0.1))
})
