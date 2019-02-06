\dontrun{
  

library(dMod)
library(ggplot2)
library(dplyr)
  
setwd(tempdir())


# Set up reactions
f <- eqnvec() %>%
  addReaction("A", "B", "k1*A", "Production of B") %>%
  addReaction("B", "C", "k2*B", "Production of C")

# Define observables and error model
observables <- eqnvec(B_obs = "B + off_B")
errors <- eqnvec(B_obs = "sqrt((sigma_rel*B_obs)^2 + sigma_abs^2)")

# Generate dMod objects
model <- odemodel(f, modelname = "errtest", compile = FALSE, solver = "deSolve")
x     <- Xs(model, optionsSens = list(method = "lsoda"), optionsOde = list(method = "lsodes"))
g     <- Y(observables, x, 
           compile = FALSE, modelname = "obsfn")
e     <- Y(errors, g, attach.input = FALSE,
           compile = FALSE, modelname = "errfn")

# Generate parameter transformation
innerpars <- getParameters(model, g, e)
covariates <- data.frame(Aini = 1:2, row.names = c("C1", "C2"))

p <- 
  eqnvec() %>%
  define("x~x", x = innerpars) %>%
  define("x~0", x = c("B", "C")) %>%
  branch(table = covariates) %>%
  insert("A~Aini", Aini = Aini) %>%
  insert("x~exp(x)", x = innerpars) %>%
  P(modelname = "parfn", compile = FALSE)

compile(g, x, e, p, output = "errtest_total")
#compile(g, x, e, p, cores = 4)


## Simulate data
ptrue <- c(k1 = -2, k2 = -3, off_B = -3, sigma_rel = log(.1), sigma_abs = log(.1))
times <- seq(0, 50, 1)
prediction <- (g*x*p)(times, ptrue, deriv = TRUE)
datasheet <- subset(as.data.frame(prediction, errfn = e), name == "B_obs")
datasheet$value <- datasheet$value + rnorm(length(datasheet$value), sd = datasheet$sigma)
data <- as.datalist(datasheet)

## Fit data with error model
obj <- normL2(data, g*x*p, e)
myfit <- trust(obj, ptrue, rinit = 1, rmax = 10, printIter = TRUE)
fits <- mstrust(obj, center = ptrue, sd = 3, fits = 2, cores = 2, printIter = TRUE)

mypars <- myfit$argument[-1]
myfixed <- myfit$argument[1]

profiles <- profile(obj + constraintL2(myfit$argument, 10), 
                    mypars, names(mypars), 
                    limits = c(-5, 5), 
                    fixed = myfixed, 
                    cores = length(mypars))
plotProfile(profiles)


## Compute prediction profile
datapoint <- datapointL2(name = "A", time = 10, value = "d1", sigma = .05, condition = "C1")
par <- trust(normL2(data, g*x*p, e) + datapoint, c(ptrue, d1 = 0), rinit = 1, rmax = 10)$argument

profile_pred <- profile(normL2(data, g*x*p, e) + datapoint, par, "d1", limits = c(-10, 10), stepControl = list(stop = "data"))

plot(profile_pred$prediction, profile_pred$data, type = "b")

}
