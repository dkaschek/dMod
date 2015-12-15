## Useful libraries --------------------------------------------

library(deSolve)
library(parallel)
library(dMod)

## Model Definition ------------------------------------------------------

# Read in model csv and generate equation list
reactionlist <- as.eqnlist(read.csv("topology.csv") )

# Add additional reactions
reactionlist <- addReaction(reactionlist, from = "", to = "", rate = "", description = "")

# Translate data.frame into equation list
f <- as.eqnvec(reactionlist)

# Define new observables based on ODE states
observables <- eqnvec(
  # y1 = "s*x1 + off"
)

# Set list of forcings
forcings <- c(
  # "u1", "u2", "u3", ...
  )


# Add observable ODEs to the original ODEs or use an observation function
# Choose one of the three options, or combine them
g <- Y(observables, f, compile = TRUE, modelname = "obsfn")


# Generate the model C files, compile them and return a list with func and extended.
model0 <- generateModel(f, forcings = forcings, modelname = "odefn")

## Parameter Transformations -------------------------------------------

# Define inner parameters (parameters occurring in the equations except forcings)
# Add names(observables) if addObservables(observables, f) is used
innerpars <- getSymbols(c(names(f), f, observables), exclude = c(forcings, "time"))

# Define additional parameter constraints, e.g. steady-state conditions
# Parameters (left-hand side) are replaced in the right-hand side of consecutive lines by resolveRecurrence() 
constraints <- resolveRecurrence(c(
  #p1 = "p2 + p3",
  #p4 = "p1*p5"   
  ))

# Build up a parameter transformation (constraints, log-transform, etc.)
# Start with the identity
trafo <- structure(innerpars, names = innerpars)
# Replace initial value parameters of the observables (if treated as states)
trafo <- replaceSymbols(names(observables), observables, trafo)
# Then employ the other parameter constraints
trafo <- replaceSymbols(names(constraints), constraints, trafo)
# Then do a log-transform of all parameters (if defined as positive numbers)
trafo <- replaceSymbols(innerpars, paste0("exp(", innerpars, ")"), trafo)


## Specify different conditions -----------------------------------------------------

conditions <- c(
  #"condition1", "condition2", ... 
); names(conditions) <- conditions

# Set condition-specific parameter transformations and generate p2p function
trafoL <- lapply(conditions, function(con) trafo)

specific <- c("")
trafoL <- lapply(conditions, function(con) {
  replaceSymbols(specific, paste(specific, con, sep = "_"), trafoL[[con]])
})

pL <- lapply(conditions, function(con) P(trafoL[[con]]))


# Set different forcings per condition
timesF <- seq(0, 100, by = 0.1)
uL <- list(
  data.frame(name = "u1", time = timesF, value = 1*dnorm(timesF, 0, 5)),
  data.frame(name = "u2", time = timesF, value = 3*dnorm(timesF, 0, 5)),
  data.frame(name = "u3", time = timesF, value = 8*dnorm(timesF, 0, 5))
); names(uL) <- conditions


# Specify prediction functions for the different conditions (depends on different forces 
# but not on different parameter transformations)

xL <- lapply(conditions, function(con) Xs(model0, uL[[con]]))

# Function for the total model prediction, returns a list of predictions (choose one of
# the two possibilities)
x <- prdfn({
  pinner <- pL[[condition]](pars, fixed)
  prediction <- xL[[condition]](times, pinner, ...)
  observation <- g(prediction, pinner, attach.input = TRUE)
  return(observation)
}, conditions = conditions)

## Data ----------------------------------------------------------------------

datasheet <- read.table("datafile.csv") # with columns name, time, value, sigma, condition
data <- as.datalist(datasheet)

## Objective Functions -------------------------------------------------------

# Initalize parameters 
outerpars <- getSymbols(do.call(c, trafoL[conditions]))
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
pouter <- rnorm(length(prior), prior, 1); names(pouter) <- outerpars

# Objective function for trust()
obj <- objfn({
  conStraintL2(pouter, prior, sigma = 10, fixed = fixed)
}, data = data, x = x, pouter = pouter, conditions = conditions)


## Howto proceed -------------------------------------------------

# Predicting and plotting
timesD <- sort(unique(unlist(lapply(data, function(d) d$time))))
times <- seq(min(timesD), max(timesD), len = 100)
prediction <- x(times, pouter)
plot(prediction)
plot(prediction, NULL, name %in% names(observables))

# Fitting
plot(data)
myfit <- trust(obj, pouter, rinit = 1, rmax = 10, iterlim = 500)
prediction <- x(times, myfit$argument)
plot(prediction, data)
plot(prediction, data, name %in% names(observables))
plot(prediction, data, name %in% names(observables), facet = "grid")

# Fitting from random positions
center <- pouter
fitlist <- mstrust(obj, center, fits = 20, cores = 4)
partable <- as.parframe(fitlist)
plotValues(partable)
bestfit <- as.parvec(partable)

# Compare predictions of best fits
plotArray(partable[1:10, ], x = x, times = times, data = data)

# Profile likelihood
profiles.approx <- do.call(rbind, mclapply(names(bestfit), function(n) profile(obj, bestfit, n, method = "integrate"), mc.cores = 4))
profiles.exact  <- do.call(rbind, mclapply(names(bestfit), function(n) profile(obj, bestfit, n, method = "optimize"), mc.cores = 4))
plotProfile(profiles.approx, profiles.exact)
plotPaths(profiles.approx, whichPar = 1)


