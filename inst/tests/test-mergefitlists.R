context("Merge multiple mfitlists")

# Set up environment
setwd("runDMod/")
source("constantPool.R")
setwd("../")

mfitlist <- mstrust(obj, pinit, fixed = fixed, fits = 2)
mfitlist3 <- rbind.fitlist(mfitlist, mfitlist, mfitlist)

attrNamesFl <- sort(names(attributes(mfitlist)))
attrNamesFl3 <- sort(names(attributes(mfitlist3)))


test_that("Short argument list workl", {
  expect_that(rbind.fitlist(), equals(NULL))
  expect_that(rbind.fitlist(mfitlist), equals(mfitlist))
})


test_that("Merging of two and three mfitlist works", {
  expect_that(max(mfitlist3$index), equals(6))
  expect_that(mfitlist3[1,-1], is_equivalent_to(mfitlist3[3,-1]))
  expect_that(mfitlist3[1,-1], is_equivalent_to(mfitlist3[5,-1]))
  expect_that(mfitlist3[2,-1], is_equivalent_to(mfitlist3[4,-1]))
  expect_that(mfitlist3[2,-1], is_equivalent_to(mfitlist3[6,-1]))
  expect_that(mfitlist3, is_a("data.frame"))
})


test_that("Attributes are correct", {
  expect_that(attrNamesFl, equals(attrNamesFl3))
  expect_that(attr(mfitlist, "metanames"), equals(attr(mfitlist3, "metanames")))
  expect_that(attr(mfitlist, "parameters"), equals(attr(mfitlist3, "parameters")))
})


fl <- attr(mfitlist, "fitlist")
fl3 <- attr(mfitlist3, "fitlist")
test_that("Attribute mfitlist is correct", {
  expect_that(length(fl3), equals(6))
  expect_that(fl[[1]]$gradient, equals(fl3[[3]]$gradient))
  expect_that(fl[[2]]$value, equals(fl3[[4]]$value))
  expect_that(fl[[1]]$argument, equals(fl3[[5]]$argument))
  expect_that(fl[[2]]$iterations, equals(fl3[[6]]$iterations))
})
  

mfitlisterrorMeta <- mfitlist
mfitlisterrorPara <- mfitlist
mfitlisterrorClass <- unclass(mfitlist)
attr(mfitlisterrorMeta, "metanames") <- c("value", "converged", "iterations")
attr(mfitlisterrorPara, "parameters") <- c("a", "b")
test_that("Errors are handled", {
  expect_that(rbind.fitlist(mfitlist, mfitlisterrorMeta), throws_error("Can not merge fitlists with differing metanames or parameters."))
  expect_that(rbind.fitlist(mfitlist, mfitlisterrorPara), throws_error("Can not merge fitlists with differing metanames or parameters."))
  expect_that(rbind.fitlist(mfitlisterrorClass), throws_error("Only data.frames can be bound"))
})