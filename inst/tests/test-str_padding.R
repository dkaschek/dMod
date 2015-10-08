context("String padding")


# Set up environment
src <- "1234567890"

paddingBlank <- paste0(rep(" ", 10), collapse = "")
paddingDots <- paste0(rep(".", 10), collapse = "")


# Access not exported function
strpad <- function(...) {strpad(...)}
environment(strpad) <- asNamespace('dMod')


test_that("string padding with default filling works properly", {
  expect_that(strpad(src, 20), equals(paste0(src, paddingBlank)))
  expect_that(strpad(src, 20, where = "right"), equals(paste0(src, paddingBlank)))
  expect_that(strpad(src, 20, where = "left"), equals(paste0(paddingBlank, src)))
})
  
test_that("string padding with custom filling works properly", {
  expect_that(strpad(src, 20, padding = '.'), equals(paste0(src, paddingDots)))
  expect_that(strpad(src, 20, padding = '.', where = "right"), equals(paste0(src, paddingDots)))
  expect_that(strpad(src, 20, padding = '.', where = "left"), equals(paste0(paddingDots, src)))
})

test_that("padding with auto elide works properly", {
  expect_that(strpad(src, 6, autoelide = TRUE), equals("123..."))
  expect_that(strpad(src, 6, where = "right", autoelide = TRUE), equals("123..."))
  expect_that(strpad(src, 6, where = "left", autoelide = TRUE), equals("...890"))
})

test_that("error are reported properly", {
  src <- "Wrong class"
  class(src) <- "MyClass"
  expect_that(strpad(src, 20), gives_warning())
  expect_that(strpad(src, 20), equals(NULL))
              
  expect_that(strpad("1", padding = ".."), throws_error())
})