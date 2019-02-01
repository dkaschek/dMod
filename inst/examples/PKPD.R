
library(dMod)
library(dplyr)
setwd(tempdir())


## Basics ------------------------------------------------------------------

# ODEs can be defined as character vectors
ODEs <- c(
  Ad = "-ka*Ad + Fabs1*INPUT1",
  Ac = "ka*Ad - Q1/Vc*Ac + Q1/Vp1*Ap1 - CL/Vc*Ac",
  Ap1 = "Q1/Vc*Ac - Q1/Vp1*Ap1",
  PL = "GR - EMAX*(Ac/Vc+eps)^hill/((Ac/Vc+eps)^hill+EC50^hill)",
  INPUT1 = "0"
)

# Dosing is encoded as events on the state INPUT1
dosing <- eventlist() %>% 
  addEvent("INPUT1", 0, "RATE", "replace") %>% 
  addEvent("INPUT1", "TINF", 0, "replace")

# ODEs and events are turned into a compiled model
# The argument "estimate" should be used to save computational time
# when the model contains many model parameters but only a few of them
# are estimated. Otherwise, sensitivity equations will be computed for all
# possible derivatives.
model <- odemodel(ODEs, events = dosing, modelname = "PKPD", 
                  estimate = c("EC50", "EMAX", "GR", "hill", "PL", "RATE"))

# A prediction function with sensitivities can be built from the model
prdfn <- Xs(model)


# The prediction function needs to be called with times and parameters
times <- seq(0, 500, 1)

inits <- c(
  Ad = 0, 
  Ac = 0, 
  Ap1 = 0, 
  PL = 10, 
  INPUT1 = 0
)

pars_PK <- c(
  Fabs1 = 1,
  ka = 0.2, 
  CL = 60,  
  Vc = 150, 
  Q1 = 10,  
  Vp1 = 790
)

pars_PD <- c(
  EC50 = 0.005,
  EMAX = 0.2,
  GR = 0.08,
  hill = 2.5
)

pars_aux <- c(eps = 1e-6)

pars_DOSING <- c(
  RATE = 600/0.1,
  TINF = 0.1
)

# Compute the prediction
prediction <- prdfn(times, c(inits, pars_PK, pars_PD, pars_DOSING, pars_aux)) 

# Plot the prediction
plot(prediction)

# Get the sensitivities and plot all sensitivities beginning with "PL"
plot(x = getDerivs(prediction), data = NULL, grepl("^PL", name))


## Extension ------------------------------------------------------------------

# It is possible to choose any parameterization of the model parameters
# The parameterization is defined as character vector
# dMod provides functions define() and insert() to help the user
# setting up transformations

# The following reparameterization inserts values,
# reparameterizes RATE by AMT and TINF
# and performs a log-transform for some parameters
trafo <- eqnvec() %>% 
  define("x~x", x = getParameters(model)) %>% 
  insert("x~0", x = c("Ad", "Ac", "Ap1", "INPUT1")) %>% 
  insert("x~y", x = names(pars_PK), y = pars_PK) %>% 
  insert("x~y", x = names(pars_aux), y = pars_aux) %>% 
  insert("RATE~AMT/TINF") %>% 
  insert("TINF~1e-4") %>% 
  insert("x~exp(x)", x = c("EC50", "EMAX", "GR", "hill"))

# The symbolic transformation is turned into a parameter transformation function
pfn <- P(trafo)

# The new parameters are a mixture of log- and non-log parameters
pars <- c(
  PL = 10,
  AMT = 600,
  GR = log(0.08),
  EC50 = log(0.005),
  EMAX = log(0.2),
  hill = log(2.5)
)

# Prediction and parameter transformation function can be concatenated to
# yield a new predictio function via the "*" operator

prediction <- (prdfn*pfn)(times, pars)
plot(prediction)
plot(x = getDerivs(prediction), data = NULL, grepl("^PL", name))



