##
## JAK-STAT pathway
##
  
# Load required packages
  library(ggplot2)
  library(deSolve)
  library(cOde)
  library(dMod)
  library(trust)
  library(parallel)
  
## Prepare equations ----------------------------------------
  setwd("~/dMod/testing/JakStat")
# Generate equations from csv  
  eq <- generateEquations(read.csv("topology.csv"))
# Define observables  
  observables <- c(tSTAT = "STAT + pSTAT + 2*pSTATdimer + off_tSTAT",
                   tpSTAT = "s_STAT*(pSTAT + 2*pSTATdimer) + off_tpSTAT")
# Define forcings and fixed states
  forcings <- "pEpoR"
  fixedStates <- c("pSTAT", "pSTATdimer", "npSTATdimer", "nSTAT1", "nSTAT2", "nSTAT3", "nSTAT4", "nSTAT5")
# Add observables to ODE  
  eq <- addObservable(observables, eq)
  model <- generateModel(eq, einspline=FALSE, forcings = forcings, fixed = fixedStates, nGridpoints = 101, jacobian = "jacvec.lsodes")
  
## Get data --------------------------------------------------
  
# Read in data.frame and plot
  mydata <- read.csv("pnas_data_original.csv")
  plotData(list(Epo = mydata)) + expand_limits(y=0) + geom_line()
# Extract forcings (pEpoR) from data  
  mydataReceptor <- subset(mydata, name=="pEpoR")
  myforc <- mydataReceptor[,c("name", "time", "value")]
  mydata <- subset(mydata, name!="pEpoR")  
  
## Parameter transformation ------------------------------------
  
# Collect all inner parameters  
  innerpars <- getSymbols(c(eq, names(eq), observables, names(observables)), exclude=forcings)
  names(innerpars) <- innerpars
  
# Define initial values / steady state value
  steadyStates <- rep("0", length(fixedStates)); names(steadyStates) <- fixedStates
  
# Define transformations
  trafo <- replaceSymbols(names(observables), observables, innerpars) #observable initial values
  trafo <- replaceSymbols(names(steadyStates), steadyStates, trafo) # steady states
  trafo <- replaceSymbols(innerpars, paste0("exp(", innerpars, ")"), trafo) #log-transform
  
# Get outer parameters from trafo
  outerpars <- getSymbols(trafo)
  
# Generate parameter transformation function
  p <- P(trafo)
  
# Initialize outer parameters  
  pini <- rep(0, length(outerpars))
  names(pini) <- outerpars
  
## Model prediction ----------------------------------------
  
# Collect times from data and augment by additional time points
  timesD <- sort(unique(mydata$time))
  times <- seq(min(timesD), max(timesD), len=250)
# Generate model prediction function
  x <- Xs(model$func, model$extended, forcings = myforc)
# Evaluate model at initial parameters and plot
  out <- x(times, p(pini))
  plotPrediction(list(states=out))
  
  
  
## Fit data ------------------------------------------------

# Define objective function
  prior <- rep(0, length(outerpars)); names(prior) <- outerpars
  myfn <- function(pp, fixed=NULL, deriv=TRUE) 
    wrss(res(mydata, x(timesD, p(pp, fixed=fixed), deriv = deriv))) + constraintL2(c(pp, fixed), prior, 20)
  
# Fit the data by a trust region algorithm
  myfit <- trust(myfn, pini + rnorm(length(pini), 0, .1), rinit=1, rmax=10, iterlim=500)
# Evaluate model at best-fitting parameter values and plot
  prediction <- x(times, p(myfit$argument))
  plotCombined(list(states=prediction, input=long2wide(mydataReceptor)), list(states=mydata, input=mydataReceptor))
  
