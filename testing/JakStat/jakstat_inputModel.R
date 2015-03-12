##
## JAK-STAT pathway with parametrized input model
##

# Load required packages
  library(ggplot2)
  library(deSolve)
  library(R2CdeSolve)
  library(trust)
  library(parallel)

## Prepare equations ---------------------------------------

# Generate equations from csv
  eq <- generateEquations(read.csv("topology_inputModel.csv"))
# Define observables
  observables <- c(tSTAT = "STAT + pSTAT + 2*pSTATdimer + offsettSTAT",
                   tpSTAT = "sSTAT*(pSTAT + 2*pSTATdimer + offsettpSTAT)",
                   pEpoR = "y1*y2")
# Add observables to ODE
  eq <- addObservable(observables, eq)
# Augment ODE by sensitivity equations
  eqAug <- c(eq, sensitivitiesSymb(eq))
# Generate C functions
  func <- funC(eq, forcings = NULL)
  funcAug <- funC(eqAug, forcings = NULL)

## Get data ------------------------------------------------

# Read in data.frame and plot
  mydata <- read.csv("pnas_data_original.csv")
  plotData(list(Epo = mydata)) + expand_limits(y=0)

## Parameter transformation --------------------------------

# Collect all inner parameters
  innerpars <- c(attr(func, "variables"), 
                 attr(func, "parameters"), 
                 "offsettSTAT", "offsettpSTAT", "scalepEpoR")
# Define initial values / steady state value
  steadyStates <- c(
    pSTAT = "0",
    pSTATdimer = "0",
    npSTATdimer = "0",
    nSTAT1 = "0",
    nSTAT2 = "0",
    nSTAT3 = "0",
    nSTAT4 = "0",
    nSTAT5 = "0",
    y1 = "0",
    y2 = "1",
    t = "0"
  )

# Define transformations
  trafo <- innerpars; names(trafo) <- innerpars #identity
  trafo <- replaceSymbols(names(observables), observables, trafo) #observable initial values
  trafo <- replaceSymbols(names(steadyStates), steadyStates, trafo) # steady states
  trafo <- replaceSymbols(innerpars, paste0("exp(log", innerpars, ")"), trafo) #log-transform

# Get outer parameters from trafo
  outerpars <- getSymbols(trafo)

# Generate parameter transformation function
  p <- P(trafo)
# Initialize outer parameters
  logpini <- rep(0, length(outerpars))
  names(logpini) <- outerpars



## Model prediction ----------------------------------------

# Collect times from data and augment by additional time points
  times <- sort(unique(c(mydata$time, seq(0, 60, len=200))))
# Generate model prediction function
  x <- Xs(func, funcAug, myforc, optionsSens = list(method="lsodes", atol=1e-10, rtol=1e-10))
# Evaluate model at initial parameters and plot
  out <- x(times, p(logpini))
  plotPrediction(list(states=out))



## Fit data ------------------------------------------------

# Define objective function
  myfn <- function(pp, fixed=NULL, doSens=TRUE) {
    o <- wrss(res(mydata, x(times, p(pp, fixed=fixed), doSens=doSens)))
    c <- constraint(pp, names(pp), 0, 10)
    
    v <- o$value + c$value
    g <- o$gradient + c$gradient
    h <- o$hessian + c$hessian
    
    return(list(value = v, gradient = g, hessian = h))
    
  }

# Fit the data by a trust region algorithm
  myfit <- trust(myfn, logpini, rinit=0.1, rmax=0.1)
# Evaluate model at best-fitting parameter values and plot
  prediction <- x(times, p(myfit$argument))
  plotCombined(list(states=prediction), list(states=mydata))
  
# Fit the data by alternative model (without nuclear export)
  fixed <- c(logp5 = -20)
  logpiniAlt <- logpini[-which(names(logpini)%in%names(fixed))]
  myfit <- trust(function(x) myfn(x, fixed), logpiniAlt, rinit=0.1, rmax=0.1)
# Evaluate alternative model for best-fitting parameters
  predictionAlt <- x(times, p(myfit$argument, fixed))
  plotCombined(list(states=prediction, alternative=predictionAlt), list(states=mydata))


  
## Multiple fits ------------------------------------------

# Fit repeatedly for different initial parameter values
# Return data.frame with the final objective value and 
# the best-fitting parameter vector
  fitlist <- do.call(rbind, mclapply(1:20, function(i) {
    mypini <- rnorm(length(outerpars), logpini, 1)
    names(mypini) <- outerpars
    myfit <- trust(myfn, mypini, rinit=0.1, rmax=1)
    out <- data.frame(chisquare = myfit$value, as.data.frame(as.list(myfit$argument)))
    return(out)
  }, mc.preschedule=FALSE))
# Sort the data.frame by increasing objective values and plot  
  fitlist <- fitlist[order(fitlist$chisquare),]
  plotFitList(wide2long(fitlist))

  
## Sensitivity analysis ----------------------------------
  
# Sensitivities for nSTAT1
  mynames <- paste("nSTAT1", outerpars[!outerpars%in%names(fixed)], sep=".")
  out <- attr(x(times, p(myfit$argument, fixed)), "deriv")
  out <- out[,c("time", mynames)]
  plotPrediction(list(test=out))
  
