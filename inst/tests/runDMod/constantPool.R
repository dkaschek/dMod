## Based on the constant pool model reduction example designed by Daniel Kascheck

library(deSolve)
library(parallel)
library(dMod)

## Model definition-------------------------------------------------------------
# Generate reaction network
f <- eqnlist()
f <- addReaction(f, "pA", "A", "k_off * pA")
f <- addReaction(f, "A", "pA", "k_on * A * exp(-0.1*time)")
fVec <- as.eqnvec(f)

# Define new observables based on ODE states
observables <- eqnvec(
  pA_obs = "scale * pA"
)

# Generate observation function
g <- Y(observables, fVec, compile = TRUE, modelname = "obsfn")
  
# Generate the model C files
model0 <- odemodel(fVec, compile = TRUE, modelname = "odefn")


## Parameter transformations----------------------------------------------------
# Define inner parameters (parameters occurring in the equations except forcings)
innerpars <- getSymbols(c(names(fVec), as.character(fVec), observables), exclude = c("time"))

# Define additional parameter constraints, e.g. initial states
constraints <- c(
  A = "1",
  pA = "0"
)

# Build up a parameter transformation (constraints, log-transform, etc.)
# Start with the identity
trafo <- structure(innerpars, names = innerpars)
# Then employ the other parameter constraints
trafo <- replaceSymbols(names(constraints), constraints, trafo)
# Then do a log-transform of all parameters (if defined as positive numbers)
trafo <- replaceSymbols(innerpars, paste0("exp(log", innerpars, ")"), trafo)
# Get names of new parameters
outerpars <- getSymbols(trafo)

# Generate parameter transformation function
p0 <- P(trafo)


## Model prediction functions---------------------------------------------------
# Generate low-level prediction function
x0 <- Xs(model0)


# Generate higher-level prediction function
y <- prdfn({
  pinner <- p0(pars, fixed)
  prediction <- x0(times, pinner, ...)
  observation <- g(prediction, pinner, attach.input = TRUE)
}, conditions = "cond1")


## Simulate Data----------------------------------------------------------------
# Use the following parameters
pouter <- c(logk_on = log(0.01),
            logk_off = log(0.1),
            logscale = log(10))

# Constant noise level and equidistant time points
noise <- 0.1
times <- seq(0, 40, by = 1.5)
set.seed(1)

# Predict model response and add noise
prediction <- y(times, pouter)
values.pA <- prediction$cond1[, "pA_obs"]
noise.pA <- rnorm(length(values.pA), 0, noise)

data <- datalist(
  cond1 = data.frame(name = "pA_obs", 
                     time = times, 
                     value = values.pA + noise.pA, 
                     sigma = noise)
)


## Objective Function-----------------------------------------------------------
# Objective function for the optimizer

obj <- objfn(data = data, pouter = pouter, conditions = names(data), x = y)


## Prepare for fitting----------------------------------------------------------
# Initalize parameters 
fixed <- c(logscale = 5)
pinit <- pouter[setdiff(outerpars, names(fixed))]

# Fit parameters and plot prediction and data
#myfit <- trust(obj, pinit, rinit = 1, rmax = 10, fixed = fixed)
#fitlist <- mstrust(obj, pinit, fixed = fixed)
