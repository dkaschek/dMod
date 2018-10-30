context("Conflicting modelnames")
test_that("modelnames behave as expected", {
  
  # What needs to be checked
  # 1. That modelname is what goes in
  # 2. That running the same code a second time breaks the behaviour of the objects of the first run
  # 3. That compiling the same structural model into a different modelname lets both functions intact
  
  #-!Start example code
  #-! library(conveniencefunctions)
  #-! library(dMod)
  library(dplyr)
  setwd(tempdir())  
  
  ## Model definition (text-based, scripting part)
  f <- NULL %>%
    addReaction("A", "B", "k1*A", "translation") %>%
    addReaction("B",  "", "k2*B", "degradation") %>%
    as.eqnvec()
  events <- eventlist(var = "A", time = 5, value = "A_add", method = "add")
  
  x1 <- odemodel(f, events = events, modelname = "odemodel") %>% Xs
  g1 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn")
  
  conditions <- c("a", "b")
  trafo <-
    getParameters(g1,x1) %>%
    setNames(.,.) %>%
    branch(conditions = conditions) %>%
    insert("x~x_cond", x = "s1", cond = condition) %>%
    insert("x~exp(x)", x = getSymbols(mytrafo[[i]])) %>%
    {.}
  
  p1 <- P(trafo, modelname = "p", compile = T)
  
  parameters <- getParameters(p1)
  pars <- structure(rnorm(length(parameters)), names = parameters)
  (g1*x1*p1)(0:10, pars)
  
  #-!End example code
  # 2. Rerunning the same parts breaks existing objects
  g2 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn")
  p2 <- P(trafo, modelname = "p", compile = T)
  
  # 3. Compiling the same structural model into a different modelname lets both functions intact
  g3 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn3")
  x3 <- odemodel(f, events = events, modelname = "odemodel3") %>% Xs
  
  # Define your expectations here
  # 1. Modelname is what goes in
  expect_equal(modelname(x1), "odemodel")
  expect_equal(modelname(g1), "obsfn")
  expect_equal(modelname(p1), paste("p", conditions, sep = "_"))
  # 2. Rerunning the same parts breaks existing objects
  expect_true(!inherits(try((g2*x1*p2)(0:10,pars)), "try-error"))
  expect_error((g1*x1*p2)(0:10, pars))
  expect_error((g2*x1*p1)(0:10, pars))
  # 3. Compiling the same structural model into a different modelname lets both functions intact
  expect_true(!inherits(try((g2*x3*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g3*x3*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g3*x1*p2)(0:10,pars)), "try-error"))
})
