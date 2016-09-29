\dontrun{

################################################################  
## Walk through the most frequently used functions in
## connection with ODE models
################################################################  
  
library(deSolve)

## Model definition (text-based, scripting part)
r <- NULL
r <- addReaction(r, "A", "B", "k1*A", "translation")
r <- addReaction(r, "B",  "", "k2*B", "degradation")

f <- as.eqnvec(r)
observables <- eqnvec(Bobs = "s1*B")

innerpars <- getSymbols(c(names(f), f, observables))
trafo0 <- structure(innerpars, names = innerpars)
trafo0 <- replaceSymbols(innerpars, 
                         paste0("exp(", innerpars, ")"), 
                         trafo0)

conditions <- c("a", "b")

trafo <- list()
trafo$a <- replaceSymbols("s1", "s1_a", trafo0)
trafo$b <- replaceSymbols("s1", "s1_b", trafo0)

events <- list()
events$a <- data.frame(var = "A", time = 5, value = 1, method = "add")
events$b <- data.frame(var = "A", time = 5, value = 2, method = "add")


## Model definition (generate compiled objects 
## and functions from above information) 

# ODE model
model <- odemodel(f)

# Observation function
g <- Y(observables, f, compile = TRUE, modelname = "obsfn")

# Parameter transformation for steady states
pSS <- P(f, method = "implicit", condition = NULL)

# Condition-specific transformation and prediction functions
p0 <- x <- NULL
for (C in conditions) {
  p0 <- p0 + P(trafo[[C]], condition = C)
  x <- x + Xs(model, events = events[[C]], condition = C)
}

## Process data

# Data
data <- datalist(
  a = data.frame(time = c(0, 2, 7),
                 name = "Bobs",
                 value = c(.3, .3, .3),
                 sigma = c(.03, .03, .03)),
  b = data.frame(time = c(0, 2, 7),
                 name = "Bobs",
                 value = c(.1, .1, .2),
                 sigma = c(.01, .01, .02))
)

## Evaluate model/data

# Initialize times and parameters
times <- seq(0, 10, .1)
outerpars <- attr(p0, "parameters")
pouter <- structure(rnorm(length(outerpars)), 
                    names = outerpars)

# Plot prediction
plot((g*x*p0)(times, pouter))
plot((g*x*pSS*p0)(times, pouter))
plotFluxes(pouter, g*x*p0, times, getFluxes(r)$B)

# Fit data with L2 constraint
constr <- constraintL2(mu = 0*pouter, sigma = 5)
myfit <- trust(normL2(data, g*x*p0) + constr, 
               pouter, rinit = 1, rmax = 10)
plot((g*x*p0)(times, myfit$argument), data)

# List of fits
obj <- normL2(data, g*x*p0) + constr
mylist <- mstrust(normL2(data, g*x*p0) + constr, 
                  pouter, fits = 10, cores = 4, sd = 4)
plot(mylist)
pars <- as.parframe(mylist)
plotValues(pars)

bestfit <- as.parvec(pars)

# Compute profiles of the fit
profile <- profile(normL2(data, g*x*p0) + constr, 
                   bestfit, "k1", limits = c(-5, 5))
plotProfile(profile)

# Add a validation objective function
vali <- datapointL2("A", 2, "mypoint", .1, condition = "a")
# Get validation point parameters
mypoint <- trust(normL2(data, g*x*p0) + constr + vali, 
                 c(mypoint = 0), rinit = 1, rmax = 10, 
                 fixed = bestfit)$argument
# Compoute profile and plot
profile <- profile(normL2(data, g*x*p0) + constr + vali, 
                   c(bestfit, mypoint), "mypoint")
plotProfile(profile)


## Using the controls() function to manipulate objects ------------------------

# List the available controls of an object
controls(x)

# List a specific control
controls(x, condition = "a", name = "optionsSens")

# Change specific control
controls(x, condition = "a", name = "optionsSens") <- 
  list(method = "lsoda", rtol = 1e-8, atol = 1e-8)

# Almost every function (objfn, parfn, prdfn, obsfn) saves the arguments 
# by which it was invoked in a list named "controls", which can be manipulated
# by the controls() function.


## New example: condition-specific observation parameters ---------------------

# Condition-specific observation parameters
f <- eqnvec(A = "-k1*A", B = "k1*A - k2*B")
observables <- eqnvec(Bobs = "s1*B")
conditions <- c("scale1", "scale2")

dynpars <- getSymbols(c(names(f), f))
obspars <- setdiff(getSymbols(observables), dynpars)
trafo0 <- structure(c(dynpars, obspars), names = c(dynpars, obspars))

obspars_all <- outer(obspars, conditions, paste, sep = "_")
trafo_all <- structure(c(dynpars, obspars_all), 
                       names = c(dynpars, obspars_all))
trafo_all <- replaceSymbols(dynpars, 
                            paste0("exp(", dynpars, ")"), 
                            trafo_all)
trafo_all <- replaceSymbols(obspars_all, 
                            paste0("exp(", obspars_all, ")"), 
                            trafo_all)

pObs <- NULL
for (C in conditions) {
  pObs <- pObs + P(replaceSymbols(obspars, 
                                  paste(obspars, C, sep = "_"), 
                                  trafo0), 
                   condition = C)
}


x <- Xs(model)
p <- P(trafo_all)
y <- g*pObs*x*p

times <- seq(0, 10, .1)
outerpars <- attr(p, "parameters")
pouter <- structure(rnorm(length(outerpars)), names = outerpars)

# Plot prediction
plot((g*pObs*x*p)(times, pouter))


}
