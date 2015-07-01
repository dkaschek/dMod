context("String eliding")


# Set up environment
src <- "1234567890"


test_that("basic eliding works for all positions", {
  expect_that(strelide(src, 6), equals("123..."))
  expect_that(strelide(src, 6, where = "right"), equals("123..."))
  expect_that(strelide(src, 6, where = "middle"), equals("12...0"))
  expect_that(strelide(src, 6, where = "left"), equals("...890"))
})

test_that("short right eliding works", {
  expect_that(strelide(src, 4, where = "right"), equals("1..."))
  expect_that(strelide(src, 3, where = "right"), equals("1.."))
  expect_that(strelide(src, 2, where = "right"), equals("1."))
  expect_that(strelide(src, 1, where = "right"), equals("."))
})

test_that("short middle eliding works", {
  expect_that(strelide(src, 5, where = "middle"), equals("1...0"))
  expect_that(strelide(src, 4, where = "middle"), equals("1..."))
  expect_that(strelide(src, 3, where = "middle"), equals("1.."))
  expect_that(strelide(src, 2, where = "middle"), equals("1."))
  expect_that(strelide(src, 1, where = "middle"), equals("."))
})

test_that("short left eliding works", {
  expect_that(strelide(src, 4, where = "left"), equals("...0"))
  expect_that(strelide(src, 3, where = "left"), equals("..0"))
  expect_that(strelide(src, 2, where = "left"), equals(".0"))
  expect_that(strelide(src, 1, where = "left"), equals("."))
})

test_that("force eliding works", {
  expect_that(strelide(src, 10, where = "right", force = TRUE), equals("1234567...")) 
  expect_that(strelide(src, 10, where = "middle", force = TRUE), equals("1234...890")) 
  expect_that(strelide(src, 10, where = "left", force = TRUE), equals("...4567890")) 
})
  
test_that("errors are reported properly", {
  src <- "Wrong class"
  class(src) <- "MyClass"
  expect_that(strelide(src, 3), gives_warning())
  expect_that(strelide(src, 3), equals(NULL))
})