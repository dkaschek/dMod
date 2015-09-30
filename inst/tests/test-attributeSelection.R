context("Attribute selection")

# Set up environment
fitlist <- readRDS(file = "dataSets/msrestore/reference.RDA")
myNames <- c("names", "row.names", "class")

test_that("Correct attributes are striped from a fitlist", {
  expect_that(attrs(fitlist), is_equivalent_to(fitlist))
  expect_that(attrs(fitlist, keep = TRUE), is_equivalent_to(fitlist))
  expect_that(attrs(fitlist, atr = "fitlist", keep = FALSE), is_equivalent_to(fitlist))
  expect_that(names(attributes(attrs(fitlist))), equals(myNames))
  expect_that(names(attributes(attrs(fitlist, atr = "fitlist", keep = FALSE))), equals(myNames))
})