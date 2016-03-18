test_that("unit test files are correctly parsed", {

  exmpextr("exampleParser-1-Input.R", "etc", "etc")
  parserOutput1 <- readLines("etc/leParser-1-Input.R")
  parserTarget1 <- readLines("etc/exampleParser-1-Expectation.R")
  
  exmpextr("exampleParser-2-Input.R", "etc", "etc")
  parserOutput2 <- readLines("etc/leParser-2-Input.R")
  parserTarget2 <- readLines("etc/exampleParser-2-Expectation.R")
  
  unlink("etc/leParser-1-Input.R")
  unlink("etc/leParser-2-Input.R")

  expect_equal(parserOutput1, parserTarget1)
  expect_equal(parserOutput2, parserTarget2)
  })
