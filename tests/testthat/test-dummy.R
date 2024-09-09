library(testthat)

test_that("start with zero, default increment", {
  input     <- 0L
  expected  <- 1L
  actual    <- dummy(input)
  expect_equal(actual, expected)
})

test_that("start with zero, default increment", {
  input     <- 0L
  increment <- 3L
  expected  <- 3L
  actual    <- dummy(input, increment)
  expect_equal(actual, expected)
})


test_that("passing a non-integer", {
  input <- "letters"
  expected_error_message <- "Must be of type 'integer'"

  expect_error(
    object    = dummy(input),
    regexp    = expected_error_message
  )
})
