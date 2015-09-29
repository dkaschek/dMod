context("Strip extra attributes from a fitlist")

# Set up environment
fitlist <- readRDS(file = "dataSets/msrestore/reference.RDA")
stripedList <- stripfitlist(fitlist)
myNames <- c("names", "row.names", "class")

test_that("original and restored fitlist are equal", {
  expect_that(stripedList, is_equivalent_to(fitlist))
  expect_that(names(attributes(stripedList)), equals(myNames))
})
