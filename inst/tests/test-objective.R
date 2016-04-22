test_that("objective functions can be evaluated and summed", {
  
  #-!Start example code
  ## Generate three objective functions
  prior <- structure(rep(0, 5), names = letters[1:5])
  
  obj1 <- constraintL2(mu = prior, attr.name = "center")
  obj2 <- constraintL2(mu = prior + 1, attr.name = "right")
  obj3 <- constraintL2(mu = prior - 1, attr.name = "left")
  
  ## Evaluate first objective function on a random vector
  pouter <- prior + rnorm(length(prior))
  #-!print(obj1(pouter))
  
  ## Split into fixed and non-fixed part
  fixed <- pouter[4:5]
  pouter <- pouter[1:3]
  #-!print(obj1(pouter, fixed = fixed))
  
  
  ## Visualize the result by a parameter profile
  myfit <- trust(obj1, pouter, rinit = 1, rmax = 10, fixed = fixed)
  myprof <- profile(obj1, myfit$argument, "a", fixed = fixed)
  #-!plotProfile(myprof)
  
  
  ## Create new objective function by adding the single ones,
  ## then evalue the random vector again
  pouter <- prior + rnorm(length(prior))
  obj <- obj1 + obj2 + obj3
  #-!print(obj(pouter))
  #-!End example code
  
  
  expect_equal(obj1(prior)$value, 0)
  expect_equal(obj2(prior)$value, 5)
  expect_equal(obj3(prior)$value, 5)
  
  expect_true(myfit$converged)
  
  expect_false(is.nan(obj(pouter)$value))
  })
