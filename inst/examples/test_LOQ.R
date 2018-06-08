library(dMod)
library(dplyr)

# Define Model and error model
x <- eqnvec(A = "-k*A") %>% odemodel(modelname = "testBLOQ") %>% Xs()
e <- eqnvec(A = "sigma_rel * A") %>% Y(x, attach.input = FALSE, compile = TRUE)

innerpars <- getParameters(x, e)
times <- seq(0, 15, .1)
p <- eqnvec() %>%
  define("pars~pars", pars = innerpars) %>% 
  define("pars~exp(pars)", pars = innerpars) %>% 
  P(condition = "C1")

# Simulate data
timesD <- c(0, 1, 2, 5, 10, 15)
parsD <- c(A = 0, k = -1, sigma_rel = 0)
data <- (x*p)(timesD, parsD) %>% 
  as.data.frame() %>% 
  mutate(sigma = pmax(0.1*value, 0.1), 
         value = rnorm(length(value), value, sigma)) %>% 
  as.datalist()


## Test for fixed sigma ----

# Test case 1: no value BLOQ
obj <- normL2(data, x*p, loq = -1)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)

# Test case 2: all values BLOQ
obj <- normL2(data, x*p, loq = 2)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)

# Test case 3: some values BLOQ
obj <- normL2(data, x*p, loq = 0.5)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)

## Test for variable sigma ----

# Test case 1: no value BLOQ
obj <- normL2(data, x*p, e, loq = -1)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)

# Test case 2: all values BLOQ
obj <- normL2(data, x*p, e, loq = 2)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)

# Test case 3: some values BLOQ
obj <- normL2(data, x*p, e, loq = 0.5)
parsT <- c(A = 0, k = -5, sigma_rel = 0)
obj(parsT)$grad
numDeriv::grad(function(x) obj(x)$value, parsT)




# Fit and plot
myfit <- trust(obj, parsT, rinit = 1, rmax = 10)
(x*p)(times, myfit$argument) %>% plot(data)
