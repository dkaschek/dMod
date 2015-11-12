context("Optimizer -- mstrust")

test_that("mstrust returns properly", {
  source("runDMod/constantPool.R")
  
  expect_is(mstrust(objfun = obj, center = pinit, studyname = "testthatmstrust", fixed = fixed), "parlist")
  
  # Clean up
  unlink(x = "testthatmstrust", recursive = TRUE)
  source("runDMod/constantPoolCleanup.R")
})
