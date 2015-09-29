context("Restore intermedite results of mstrust")

# Set up environment
fitlist <- readRDS(file = "dataSets/msrestore/reference.RDA")
restoredFitlist <- msrestore("dataSets/msrestore/interres/")

test_that("original and restored fitlist are equal", {
  expect_that(restoredFitlist, is_equivalent_to(fitlist))
})
