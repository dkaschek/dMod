context("funC0")
test_that("funC0 behaviour of compiled and uncompiled match", {
  
  #-!Start example code
  #-! library(dMod)
  
  myfun1 <- funC0(c(y = "a*x^4 + b*x^2 + c"), modelname = "fun1", compile = T)
  myfun2 <- funC0(c(y = "a*x^4 + b*x^2 + c"))
  out1 <- myfun1(a = -1, b = 2, c = 3, x = seq(-2, 2, .1), attach.input = TRUE)
  out2 <- myfun2(a = -1, b = 2, c = 3, x = seq(-2, 2, .1), attach.input = TRUE)
  
  #-!End example code
  
  
  # Define your expectations here
  expect_true(file.exists("fun1.c"))
  expect_identical(out1, out2)
  
  # remove c files
  unlink("fun1.c")
  unlink("fun1.o")
  unlink("fun1.so")
})
