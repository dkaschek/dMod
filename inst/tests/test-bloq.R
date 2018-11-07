
context("wrss/nll")
test_that("safe_numerics_for_bloq", {

  #-!Start example code
  #-! library(dMod)
  #-! library(dplyr)
  library(dplyr)
  x <- Xt()
  g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"))
  p <- P(c("a" = "a", "b" = "b"), condition = "C1")
  times <- seq(-5, 5, by = .05)
  pars <- c(a = 10, b = 1)
  sigma <- 1
  
  objvals <- list()
  
  #-!plot((g*x*p)(times, pars))
  
  set.seed(1)
  data1 <- data.frame(name = "y", time = times, sigma = sigma, condition = "C1", stringsAsFactors = F) %>% 
    mutate(value = pars[1]*time^2 + pars[2] + rnorm(length(time), sd = sigma)) %>% 
    as.data.frame() %>% as.datalist
  
  obj1 <- normL2(data1, (g*x*p))
  objvals[[1]] <- obj1(c(a = 0, b  =100))
  
  fit <- trust(obj1, pars, 1, 10)
  #-!plot((g*x*p)(times, fit$argument), data1)
  
  set.seed(1)
  data2 <- data.frame(name = "y", time = times, sigma = sigma, condition = "C1", stringsAsFactors = F) %>% 
    mutate(value = pars[1]*time^2 + pars[2] + rnorm(length(time), sd = sigma),
           lloq = 50) %>% 
    as.data.frame() %>% as.datalist
  
  obj2 <- normL2(data2, (g*x*p))
  objvals[[2]] <- obj2(c(a = 0, b  =100))
  
  
  fit2 <- trust(obj2, pars, 1, 10)
  #-!plot((g*x*p)(times, fit2$argument), data2)
  
  e <- Y(c(y = "s0"), g)
  p <- getParameters(g,x,e) %>% setNames(.,.) %>% P(condition = "C1")
  obj3 <- normL2(data1, (g*x*p), e)
  objvals[[3]] <- obj3(c(a = 0, b  =100, s0 = 0.1))
  
  obj4 <- normL2(data2, (g*x*p), e)
  objvals[[4]] <- obj4(c(a = 0, b  =100, s0 = 0.1))
  
  #-! print(objvals)
  #-!End example code
  values_finite <- sapply(objvals, function(.x) is.finite(.x$value))
  gradients_finite <- do.call(c, lapply(objvals, function(.x) is.finite(.x$gradient)))
  
  # Define your expectations here
  expect_true(all(values_finite))
  expect_true(all(gradients_finite))
})
  
