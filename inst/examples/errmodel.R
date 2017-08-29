\dontrun{
  

library(dMod)
library(ggplot2)
setwd("/tmp")

# Set up reactions
f <- NULL
f <- addReaction(f, "A", "B", "k1*A", "Production of B")
f <- addReaction(f, "B", "C", "k2*B", "Production of C")

# Define observables and error model
observables <- eqnvec(B_obs = "B + off_B")
errors <- eqnvec(B_obs = "sqrt((sigma_rel*B_obs)^2 + sigma_abs^2)")

# Generate dMod objects
model <- odemodel(f, modelname = "errtest", compile = FALSE, solver = "Sundials")
x     <- Xs(model, optionsSens = list(method = "bdf"), optionsOde = list(method = "bdf"))
g     <- Y(observables, x, 
           compile = FALSE, modelname = "obsfn")
e     <- Y(errors, g, attach.input = FALSE,
           compile = FALSE, modelname = "errfn")

# Generate parameter transformation
innerpars <- getParameters(model, g, e)
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~1", x = "A", trafo)
trafo <- repar("x~0", x = c("B", "C"), trafo)
trafo <- repar("x~exp(x)", x = innerpars, trafo)

p <- P(trafo, condition = "C1", modelname = "parfn", compile = FALSE) +
     P(trafo, condition = "C2", modelname = "parfn", compile = FALSE)

compile(g, x, e, p, output = "errtest_total")
compile(g, x, e, p, cores = 4)


## Simulate data
ptrue <- c(k1 = -2, k2 = -3, off_B = -3, sigma_rel = log(.1), sigma_abs = log(.1))
times <- seq(0, 50, 1)
prediction <- (g*x*p)(times, ptrue, deriv = TRUE)
datasheet <- subset(as.data.frame(prediction, errfn = e), name == "B_obs")
datasheet$value <- datasheet$value + rnorm(length(datasheet$value), sd = datasheet$sigma)
data <- as.datalist(datasheet)

## Fit data with error model
obj <- normL2(data, g*x*p, e)
myfit <- trust(obj, ptrue, rinit = 1, rmax = 10)
fits <- mstrust(obj, center = ptrue, sd = 3, fits = 10)
profiles <- profile(obj + constraintL2(myfit$argument, 10), 
                    myfit$argument, names(myfit$argument), limits = c(-5, 5), cores = 4)
plotProfile(profiles)

## Fit externally
out <- runbg({
  trust(obj, ptrue, rinit = 1, rmax = 10)
}, machine = "localhost", filename = "test", input = c("obj", "ptrue"), compile = TRUE)

## Fit on grid
out <- runbg_bwfor({
  trust(obj, ptrue, rinit = 1, rmax = 10)
}, machine = "bwfor", filename = "test", input = c("obj", "ptrue"), compile = TRUE, nodes = 2, cores = 1, walltime = "00:01:00")



## Plotting
out <- as.data.frame((g*x*p)(times = seq(0, 50, len = 100), pars = myfit$argument), errfn = e)
ggplot(out, aes(x = time, y = value, ymin = value-sigma, ymax = value+sigma, color = condition, fill = condition)) +
  facet_wrap(~name, scales = "free") +
  geom_line() + geom_ribbon(alpha = .2, lty = 0) + 
  geom_point(data = as.data.frame(data)) +
  theme_dMod() + scale_color_dMod() + scale_fill_dMod()


## Compute prediction profile
datapoint <- datapointL2(name = "A", time = 10, value = "d1", sigma = .05, condition = "C1")
par <- trust(normL2(data, g*x*p, e) + datapoint, c(ptrue, d1 = 0), rinit = 1, rmax = 10)$argument

profile_pred <- profile(normL2(data, g*x*p, e) + datapoint, par, "d1", limits = c(-10, 10), stepControl = list(stop = "data"))

plot(profile_pred$prediction, profile_pred$data, type = "b")

}
