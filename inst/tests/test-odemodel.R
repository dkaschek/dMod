context("odemodel")
test_that("Basic ODE model integration with forcings", {
  
  ## Generate the same model from an equation list
  f <- addReaction(NULL, from = "", to = "A", rate = "switch*F", description = "production")
  f <- addReaction(f   , from = "A", to = "", rate = "k*A", description = "degradation")
  model <- odemodel(f, forcings = "F", fixed = "switch")
  
  # create forcings
  forc1 <- data.frame(name = "F", time = seq(0,5, 1), value = sin(seq(0,5,1)))
  forc2 <- data.frame(name = "F", time = seq(0,5, 1), value = exp(-seq(0,5,1)))
  forc3 <- data.frame(name = "F", time= 0           , value = 0.1)
  
   
  
  # Default settings
  x <- 
    Xs(model, forc1, condition = "forc1") + 
    Xs(model, forc2, condition = "forc2") + 
    Xs(model, forc3, condition = "forc3")
  
  times <-  seq(0,5, 0.001)
  pars <- setNames(runif(length(getParameters(x))), getParameters(x))
  pred <- x(times, pars)
  
  # Expect that slope changes four times
  slope1 <- diff((pred$forc1)[, "F"])
  expect_true(length(which(abs(diff(slope1)) > 1e-8)) == 4)
  
  
  # Constant interpolation
  x <- 
    Xs(model, forc1, condition = "forc1", fcontrol = list(method = "constant", rule = 2)) + 
    Xs(model, forc2, condition = "forc2", fcontrol = list(method = "constant", rule = 2)) + 
    Xs(model, forc3, condition = "forc3", fcontrol = list(method = "constant", rule = 2))
  
  pred <- x(times, pars)
  
  # Expect that slope is zero except for 4 jumps
  slope1 <- diff((pred$forc1)[, "F"])
  expect_true(length(which(abs(slope1) > 1e-8)) == 4)
  
  # Constant interpolation without sensitivities
  x <- 
    Xf(model, forc1, condition = "forc1", fcontrol = list(method = "constant", rule = 2)) + 
    Xf(model, forc2, condition = "forc2", fcontrol = list(method = "constant", rule = 2)) + 
    Xf(model, forc3, condition = "forc3", fcontrol = list(method = "constant", rule = 2))
  
  pred <- x(times, pars)
  
  # Expect that slope is zero except for 4 jumps
  slope1 <- diff((pred$forc1)[, "F"])
  expect_true(length(which(abs(slope1) > 1e-8)) == 4)
  
  
})