## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 10, fig.height = 8, warning = FALSE, message = FALSE)

## ------------------------------------------------------------------------
library(deSolve)
library(dMod)
set.seed(1)

## ------------------------------------------------------------------------
# Generate the ODE model
reactions <- as.eqnlist(read.csv("topology.csv"))
print(reactions)

## ------------------------------------------------------------------------
# Parameterize the receptor phosphorylation
receptor <- "pEpoR*((1 - exp(-time*lambda1))*exp(-time*lambda2))^3" 
reactions$rates <- replaceSymbols(
  what = "pEpoR", 
  by = receptor,
  x = reactions$rates
)

## ------------------------------------------------------------------------
# Define parameters that are not estimated
fixed.zero <- setdiff(reactions$states, "STAT") 
# Generate odemodel
model0 <- odemodel(reactions, modelname = "jak-stat", compile = TRUE, 
                   fixed = fixed.zero, jacobian = "inz.lsodes")

## ---- fig.width = 6, fig.height = 4--------------------------------------
# Generate a prediction function
x <- Xs(model0, optionsSens = list(method = "lsodes"))

# Make a prediction based on random parameter values
parameters <- getParameters(x)
pars <- structure(runif(length(parameters), 0, 1), names = parameters)
times <- seq(0, 10, len = 100)
prediction <- x(times, pars)
plot(prediction)


## ------------------------------------------------------------------------
# Define observables like total STAT, total phosphorylated STAT, etc.
observables <- eqnvec(
  tSTAT = "s_tSTAT*(STAT + pSTAT + 2*pSTATdimer) + off_tSTAT",
  tpSTAT = "s_tpSTAT*(pSTAT + 2*pSTATdimer) + off_tpSTAT",
  pEpoR = paste0("s_EpoR *", receptor)
)

# Define the observation function. Information about states and dynamic parameters
# is contained in reactions
g <- Y(observables, reactions, modelname = "obsfn", compile = TRUE, attach.input = FALSE)

## ---- fig.width = 6, fig.height = 2.5------------------------------------
# Make a prediction of the observables based on random parameter values
parameters <- union(getParameters(x), getParameters(g))
pars <- structure(runif(length(parameters), 0, 1), names = parameters)
times <- seq(0, 10, len = 100)
prediction <- (g*x)(times, pars)
plot(prediction)

## ------------------------------------------------------------------------
# Start with the identity transformation
innerpars <- union(attr(x, "parameters"), attr(g, "parameters"))
trafo <- as.eqnvec(innerpars, names = innerpars)

# Fix some initial values
trafo[fixed.zero] <- "0"

# Log-transform
trafo <- replaceSymbols(innerpars, paste0("exp(", innerpars, ")"), trafo)

# Generate the parameter transformation function
p <- P(trafo, condition = "Epo")


## ------------------------------------------------------------------------
# Modify rate p1 
trafo["pEpoR"] <- paste("multiple *", trafo["pEpoR"])

# Add new parameter transformation function to the existing one
p <- p + P(trafo, condition = "Epo prediction")

## ---- fig.width = 6, fig.height = 2.5------------------------------------
# Make a prediction of the observables based on random parameter values
parameters <- getParameters(p)
pars <- structure(runif(length(parameters), 0, 1), names = parameters)
pars["multiple"] <- 2
times <- seq(0, 10, len = 100)
prediction <- (g*x*p)(times, pars)
plot(prediction)

## ---- fig.width = 6, fig.height = 2--------------------------------------
datasheet <- read.csv("pnas_data_original.csv")
data <- as.datalist(datasheet, split.by = "condition")
plot(data)

## ------------------------------------------------------------------------
obj <- normL2(data, g*x*p)

## ------------------------------------------------------------------------
mu <- structure(rep(-1, length(parameters)), names = parameters)
mu["multiple"] <- 0
constr <- constraintL2(mu = mu, sigma = 5)

## ---- fig.width = 6, fig.height = 2--------------------------------------

myfit <- trust(obj + constr, mu, rinit = 1, rmax = 10)
times <- 0:60
plot((g*x*p)(times, myfit$argument), data)


## ---- fig.width = 5, fig.height = 4--------------------------------------

set.seed(1)
fitlist <- mstrust(obj + constr, center = mu, fits = 20, cores = 1, min = -1, max = 1, samplefun = "runif")
pars <- as.parframe(fitlist)
plotValues(subset(pars, converged))
plotPars(subset(pars, converged))


## ---- fig.width = 6, fig.height = 5--------------------------------------

controls(g, NULL, "attach.input") <- TRUE
plotArray(subset(pars, converged), g*x*p, 0:60, data, !grepl("prediction", name))


## ---- fig.width = 4, fig.height = 3--------------------------------------

myprofile <- profile(obj + constr, pars = myfit$argument, whichPar = "s_EpoR")
plotProfile(myprofile)


## ---- fig.width = 6, fig.height = 5--------------------------------------

plotPaths(myprofile)


## ---- fig.width = 4, fig.height = 3--------------------------------------

fixed <- c(p1 = 0, pEpoR = 0) # log values
pars <- mu[setdiff(names(mu), names(fixed))]
myfit <- trust(obj + constr, pars, rinit = 1, rmax = 10, fixed = fixed)
myprofile <- profile(obj + constr, pars = myfit$argument, whichPar = "s_EpoR", fixed = fixed)

plotProfile(myprofile)


