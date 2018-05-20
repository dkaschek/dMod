## Generate a prediction function
regfn <- c(y = "sin(a*time)")

g <- Y(regfn, parameters = "a")
x <- Xt(condition = "C1")

## Generate data
data <- datalist(
  C1 = data.frame(
    name = "y",
    time = 1:5,
    value = sin(1:5) + rnorm(5, 0, .1),
    sigma = .1
  )
)

## Initialize parameters and time 
pars <- c(a = 1)
times <- seq(0, 5, .1)

## Do many fits from random positions and store them into parlist
out <- as.parlist(lapply(1:50, function(i) {
  trust(normL2(data, g*x), pars + rnorm(length(pars), 0, 1), rinit = 1, rmax = 10)
}))

## Reduce parlist to parframe
parframe <- as.parframe(out)
plotValues(parframe)

## Get steps
getStepIndices(parframe, nsteps = 2, tol = 1)

getSteps(parframe, nsteps = 2)

plotValues(getSteps(parframe, nsteps = 2))
