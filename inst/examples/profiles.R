## Parameter transformation
trafo <- eqnvec(a = "exp(loga)", 
                b = "exp(logb)", 
                c = "exp(loga)*exp(logb)*exp(logc)")
p <- P(trafo)

## Objective function
obj1 <- constraintL2(mu = c(a = .1, b = 1, c = 10), sigma = .6)
obj2 <- constraintL2(mu = c(loga = 0, logb = 0), sigma = 10)
obj <- obj1*p + obj2

## Initialize parameters and obtain fit
pars <- c(loga = 1, logb = 1, logc = 1)
myfit <- trust(obj, pars, rinit = 1, rmax = 10)
myfit.fixed <- trust(obj, pars[-1], rinit = 1, rmax = 10, fixed = pars[1])

## Compute profiles by integration method
profiles.approx <- do.call(
  rbind, 
  lapply(1:3, function(i) {
    profile(obj, myfit$argument, whichPar = i, limits = c(-10, 10),
            method = "integrate")
  })
)

## Compute profiles by repeated optimization 
profiles.exact <- do.call(
  rbind, 
  lapply(1:3, function(i) {
    profile(obj, myfit$argument, whichPar = i, limits = c(-10, 10),
            method = "optimize")
  })
)

## Compute profiles for fit with fixed element by integration method
profiles.approx.fixed <- do.call(
  rbind, 
  lapply(1:2, function(i) {
    profile(obj, myfit.fixed$argument, whichPar = i, limits = c(-10, 10),
            method = "integrate",
            fixed = pars[1])
  })
)

## Plotting
plotProfile(profiles.approx)
plotProfile(list(profiles.approx, profiles.exact))
plotProfile(list(profiles.approx, profiles.approx.fixed))

plotPaths(profiles.approx, sort = TRUE)
plotPaths(profiles.approx, whichPar = "logc")
plotPaths(list(profiles.approx, profiles.approx.fixed), whichPar = "logc")
