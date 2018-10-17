rm(list = ls())
setwd(tempdir())


library(dMod)
library(dplyr)
library(ggplot2)

# Define reactions
reactions <- eqnlist() %>% 
  addReaction("Ad" , "Ac" , "ka*Ad"     ) %>% 
  addReaction("Ac" , "Ap1", "Q1*Ac/Vc"  ) %>% 
  addReaction("Ap1", "Ac" , "Q1*Ap1/Vp1") %>% 
  addReaction("Ac" , "0"  , "CL*Ac/Vc"  )

# Define dosing
dosing <- eventlist() %>% 
  addEvent("Ad", time = "5 + tlag", value = "dose")

# Generate prediction function x(t, p)
x <- reactions %>% 
  odemodel(events = dosing, compile = TRUE) %>% 
  Xs()

# Define parameters and simulation time
pars <- c(
  Ad = 0  , Ac = 0  , Ap1 = 0  , 
  CL = 3, ka = 0.1, Q1  = .2, Vc = 10, Vp1 = 100,
  tlag = 5, dose = 100)

times <- seq(0, 40, .1)

# Simulate model and outputs
x(times, pars) %>% plot()



x(times, pars) %>% plot() + theme(legend.position = "none")


# Generate observation function g(x)
g <- eqnvec(OUTPUT1 = "Ac/Vc") %>% 
  Y(x, compile = TRUE)


# Simulate model and outputs
x(times, pars) %>% plot()
(g*x)(times, pars) %>% plot()





g <- eqnvec(Cc = "Ac/Vc") %>% 
  Y(x, compile = TRUE, attach.input = FALSE)



e <- eqnvec(Cc = "sigma_abs + Cc * sigma_rel") %>% 
  Y(g, modelname = "err_PK", attach.input = FALSE, compile = TRUE)

p <- eqnvec() %>% 
  define("x~x", x = getParameters(g, x, e)) %>% 
  define("x~0", x = c("Ac", "Ad", "Ap1")) %>% 
  insert("x~exp(x)", x = c("CL", "ka", "Q1", "Vc", "Vp1", "sigma_abs", "sigma_rel")) %>% 
  P(modelname = "par_PK", compile = TRUE)


# Simulate
pars <- c(CL = log(3), ka = log(0.1), Q1 = log(0.2), Vc = log(10), Vp1 = log(100),
          dose = 200, tlag = 1.5,
          sigma_abs = 0.02, sigma_rel = 0.1)

times <- seq(0, 60, .1)

prd <- g*x*p
prd(times, pars) %>% plot()

# Split up parameterization

