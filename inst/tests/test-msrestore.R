context("Restore intermedite results of mstrust")

# Set up environment
fitlist <- readRDS(file = "dataSets/msrestore/reference.RDA")
restoredFitlist <- msrestore("dataSets/msrestore/interres/")

test_that("original and restored fitlist are equal", {
  expect_that(restoredFitlist, is_equivalent_to(fitlist))
})

# Cleanup
unlink("runDMod/*.c")
unlink("runDMod/*.o")
unlink("runDMod/*.so")
unlink("runDMod/mstrust.log")
unlink("runDMod/fitlist.rda")
unlink("mstrust.log")
unlink("fitlist.rda")