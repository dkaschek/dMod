context("summation of objlists")
test_that("+.objlist produces expectable output", {
  
  #-!Start example code
  # library(dMod)
  nm1 <- letters[1:3]
  nm2 <- letters[2:4]
  
  ol1 <- objlist(1, `names<-`(rep(1,3), nm1), matrix(1, 3,3, F, list(nm1, nm1)))
  ol2 <- objlist(2, `names<-`(rep(2,3), nm2), matrix(2, 3,3, F, list(nm2, nm2)))
  # same elements, same grad and hessian names
  ol1 + ol1
  # same elements, different gradient and hessian names
  ol1 + ol2
  #-!End example code
  
  
  sum_same <- ol1 + ol1
  expectation_sum_same <- structure(list(value = 2, 
                                         gradient = structure(c(2, 2, 2), .Names = c("a", "b", "c")), 
                                         hessian = structure(c(2, 2, 2, 2, 2, 2, 2, 2, 2), 
                                                             .Dim = c(3L, 3L), .Dimnames = list(c("a", "b", "c"), c("a", "b", "c")))), 
                                    .Names = c("value", "gradient", "hessian"), class = "objlist")
  
  sum_diff <- ol1 + ol2
  expectation_sum_diff <- structure(list(value = 3, gradient = structure(c(1, 3, 3, 2), .Names = c("a",
                                                                                                   "b", "c", "d")), hessian = structure(c(1, 1, 1, 0, 1, 3, 3, 2,
                                                                                                                                          1, 3, 3, 2, 0, 2, 2, 2), .Dim = c(4L, 4L), .Dimnames = list(c("a",
                                                                                                                                                                                                        "b", "c", "d"), c("a", "b", "c", "d")))), .Names = c("value",
                                                                                                                                                                                                                                                             "gradient", "hessian"), class = "objlist")
  
  # check for well behaved objlists
  ol3 <- objlist(3, c("A" = 3), 0)
  
  
  # Define your expectations here
  expect_identical(sum_same, expectation_sum_same)
  expect_identical(sum_diff, expectation_sum_diff)
  expect_error({ol1 + ol3})
})
