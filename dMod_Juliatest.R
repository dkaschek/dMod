## Load Libraries -----
rm(list=ls())
library(dMod)
library(dplyr)
library(ggplot2)
theme_set(theme_dMod())
setwd(tempdir()) # change if you want to keep the files that dMod creates internally
set.seed(5555)
## Set up the ODE model -----------------------------------------------------------------------

# Define via reactions

reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell",  rate = "import*TCA_buffer", description = "Uptake") %>%
  addReaction("TCA_cell", "TCA_buffer",  rate = "export_sinus*TCA_cell", description = "Sinusoidal export") %>%
  addReaction("TCA_cell", "TCA_cana",    rate = "export_cana*TCA_cell", description = "Canalicular export") %>%
  addReaction("TCA_cana", "TCA_buffer",  rate = "reflux*TCA_cana", description = "Reflux into the buffer")

# Define via ODE system

# f <- eqnvec(TCA_buffer = "-import * TCA_buffer + export_sinus * TCA_cell + reflux * TCA_cana",
#             TCA_cana = "export_cana * TCA_cell - reflux * TCA_cana",
#             TCA_cell = "import * TCA_buffer - export_sinus * TCA_cell - export_cana * TCA_cell")

# Translate reactions into ODE model object
mymodel <- odemodel(reactions, modelname = "bamodel", solver = "OrdinaryDiffEq")


# Generate trajectories for the default condition
x <- Xs(mymodel, condition = NULL)


# For demonstration define parameters (initials and dynamic parameters)
pars <- c(reflux = 0.1, TCA_buffer = 1, TCA_cell = 0, TCA_cana = 0, import = 0.2, export_sinus = 0.2, export_cana = 0.04)



# Plot trjectories
times <- seq(0, 50, 0.01) 
out <- x(times, pars)
plot(out) # Note this is a ggplot object

# Define observables buffer and cellular
observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL, compile = TRUE, modelname = "bamodel")



## Simulate data ---------------------------------------------------------------------------

# Fix parameters to steady state values (calculated on paper / by hand), replace buffer -> TCA_buffer = 0
pars["TCA_cell"] <- 0.3846154
pars["TCA_cana"] <- 0.1538462
pars["TCA_buffer"] <- 0
pars["s"] <- 1e3

# Simulate trajectories 
# Evaluate the expression separately
out <- (g * x)(times, pars, conditions = "closed")
plot(out)





# simulate data for time points
timesD <- c(0.1, 1, 3, 7, 11, 15, 20, 41)
datasheet <- subset(as.data.frame(out), time %in% timesD & name %in% names(observables))
datasheet <- within(datasheet, {
  sigma <- 0.1*value
  value <- rnorm(length(value), value, sigma)
})
data <- as.datalist(datasheet)
plot(out, data)


# Define parameter transformations explicitly without convinience functions insert(), define() and branch() or repar() (old)

trafo <- eqnvec(TCA_buffer = "0", 
                TCA_cell = "10^(TCA_CELL)",
                TCA_cana = "10^(TCA_CANA)")\

pars <- c(TCA_CELL = 0.1, TCA_CANA = 0.1)

p <- P(trafo, condition = "closed")

p(pars)
getParameters(p)

# Set every parameter to -1 in the log-space
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*p)(times, pouter),data)




## Use simulate data to calibrate outer model parameters -------------------------------------------------------------

# # One Fit
obj <- normL2(data, g * x * p) + constraintL2(pouter, sigma = 20)


# Fit on time (starting from pouter)

myfit_1 <- trust(obj, pouter, rinit = 0.1, rmax = 5)
mypred_1 <- (g * x * p)(times, myfit_1$argument)
plot(mypred_1, data)

## Handling different experimental conditions

# Define reflux and open condition according to "Dynamic Modelling in R p. 19-20"
pars["reflux"] <- 1e3
out <- (g * x)(times, pars, conditions = "open")
datasheet <- subset(as.data.frame(out), time %in% timesD & name %in% names(observables))
datasheet <- within(datasheet, {
  sigma <- 0.1*value
  value <- rnorm(length(value), value, sigma)
})
data <- data + as.datalist(datasheet)
plotData(data)

# Parameter Trafo, usage of "+" operator for trafo functions (output of P())
trafo <- getEquations(p, conditions = "closed")
trafo["reflux"] <- "10^(reflux_open)"
p <- p + P(trafo, condition = "open")

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*p)(times, pouter),data)


# Objective function
obj <- normL2(data, g * x * p) + constraintL2(pouter, sigma = 20)

# Evaluation of obj at pouter
obj(pouter)


# Fit 50 times, sample with sd=4 around pouter
out_frame <- mstrust(obj,pouter, sd = 4,  studyname = "bamodel", cores=4,fits=20, iterlim = 200)
out_frame <- as.parframe(out_frame)
plotValues(out_frame) # Show "Waterfall" plot
plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(out_frame)


# Plot predictions along data
plot((g * x * p)(times, bestfit), data)

# Plot sensis
plot(getDerivs((g * x * p)(times, bestfit)))

# Calculate Parameter Profiles and plot different contributions (for identifiablility only "data" is of interest)
profiles <- profile(obj, bestfit, names(bestfit), method = "optimize", cores = 12)
plotProfile(profiles)
plotProfile(profiles,mode %in% c("data", "prior"))
plotPaths(profiles, whichPar = "TCA_cana")
plotPaths(profiles, whichPar = "export_cana")
plotPaths(profiles, whichPar = "s")

# Tighten model assumptions with steady state constraint
pSS <- NULL
equations <- getEquations(p) # getEquations retrieves the "trafo" from p 
conditions <- names(equations)

# Set up steady state constraint
for (n in conditions) {
  equations[[n]]["TCA_cana"] <- "10^(export_cana)*10^(TCA_cell)/10^(reflux)"
  pSS <- pSS + P(equations[[n]], condition = n)
}
outerpars <- getParameters(pSS)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*pSS)(times, pouter),data)

# Objective function
obj <- normL2(data, g * x * pSS, times = times, attr.name = "data") + constraintL2(pouter, sigma = 5, attr.name = "prior")

# Fit 50 times again
out_frame <- mstrust(obj,pouter,studyname = "bamodel", cores=10,fits=50)
out_frame <- as.parframe(out_frame)
plotValues(out_frame) # Show "Waterfall" plot
plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(out_frame)
plot((g * x * pSS)(times, bestfit), data)


# Calculate Parameter Profiles
profiles <- profile(obj, bestfit, names(bestfit), cores = 10, limits = c(lower = -10, upper = 10),
                    stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 100, stop = "data"))
plotProfile(profiles,mode %in% c("data", "prior"))
plotPaths(profiles, whichPar = "s")


# s and initial value TCA_cell render model structurally non identifiable
#  --> fix s to 1 (on log scale to 0)

trafo <- getEquations(pSS)
trafo <- repar("x~0", x = "s", trafo)
p <- P(trafo)


# Objective function
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
obj <- normL2(data, g * x * p, times = times, attr.name = "data") + constraintL2(pouter, sigma = 20, attr.name = "prior")

# Fit 50 times again
out_frame <- mstrust(obj,pouter,studyname = "bamodel", cores=10,fits=50)
out_frame <- as.parframe(out_frame)
plotValues(out_frame) # Show "Waterfall" plot
plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(out_frame)
plot((g * x * p)(times, bestfit), data)

# remove the prior
obj <- normL2(data, g * x * p, times = times, attr.name = "data")

# refit
refit <- trust(obj, bestfit, rinit = 1, rmax = 10)
plot((g * x * p)(times, refit$argument), data)

profiles <- profile(obj, refit$argument, names(refit$argument), cores = 10, limits = c(lower = -5, upper = 5), method = "optimize",
                    stepControl = list(stop = "data"))
plotProfile(profiles, mode == "data")
plotPaths(profiles, whichPar = "reflux_open")

# Note: The parameter "reflux_open" is still practical non-identifiable 
# The profile is open to the left -> A possible reduction of the model would be the assumption of an immediate reflux, i.e. the limit case reflux_open -> infinity
# An pragmatic but ugly way to circumvent the reformulation of the ODE system is to fix the reflux_open parameter to a high value. E.g. 1e3

# One could also check the models ability to produce reliable predictions, by the calculation of prediction uncertainty with profile likelihood
# The calculation of prediction confidence intervals is done at next

## Prediction uncertainty taken from validation profile --------------------------------------------------------------------------

# choose sigma below 1 percent of the prediction in order to pull the prediction strongly towards d1
obj.validation <- normL2(data, g * x * p, times = times, attr.name = "data") +
  datapointL2(name = "TCA_cell", time = 10, value = "v", sigma = 1, attr.name = "validation", condition = "closed")

# If sigma is not known, and you therefore decide to calculate prediction confidence intervals, just choose a very small sigma, in order to "pull strongly" on the trajectory

# refit
myfit <- trust(obj.validation, parinit = c(v = 180, bestfit), rinit = 1, rmax = 10, iterlim = 1000)

# Calculate profile
validation_profile <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "optimize",
                              optControl = list(rinit = .1, rmax = 10, iterlim = 100))


# plotProfile(validation_profile) # This also plos the prediction colums, which is a bug in the code.
plotProfile(validation_profile, mode %in% c("validation", "data")) # Plots only the two contributions validation and data, along with the sum (total)
# Is the contribution of validation small?? # If yes: The total aligns with a prediction profile

# Confidence Interval of the prediction
confint(validation_profile, val.column = "value")



## Prediction band (prediction uncertainty for several time points) --------------------------------------------------------------
# Here we calculate a prediction CI for different timepoints. In the end we interpolate to a "prediction band"
prediction_band <- do.call(rbind, lapply(seq(10, 50, 10), function(t) {
  
  cat("Computing prediction profile for t =", t, "\n")
  
  obj.validation <- normL2(data, g * x * p, times = times, attr.name = "data") + 
    datapointL2(name = "TCA_cell", time = t, value = "v", sigma = 1, attr.name = "validation", condition = "closed")
  
  refit <- trust(obj.validation, parinit = c(v = 180, bestfit), rinit = 1, rmax = 10, iterlim = 1000)
  
  profile_prediction <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "optimize")
  
  d1 <- confint(profile_prediction, val.column = "value")
  
  # Output
  data.frame(time = t, condition = "closed", name = "TCA_cell",  d1[-1])
  
}))

prediction <- (g * x * p)(times, refit$argument) %>%
  as.data.frame() 


prediction_band_spline <- data.frame(
  time = prediction$time[prediction$time>=10],
  value = prediction$value[prediction$time>=10],
  condition = "closed",
  name = "TCA_cell",
  lower = spline(prediction_band$time, prediction_band$lower, xout = prediction$time[prediction$time>=10])$y,
  upper = spline(prediction_band$time, prediction_band$upper, xout = prediction$time[prediction$time>=10])$y
)

# Create the ggplot
ggplot(prediction, aes(x = time, y = value, color = condition)) +
  geom_line() +  # Line connecting the points for each condition
  geom_ribbon(data = prediction_band_spline, aes(x = time, ymin = lower, ymax = upper, fill = condition), 
              lty = 0, alpha = .3, show.legend = F) +  # Show ribbon in the legend
  geom_point(data = prediction_band, aes(x = time, y = lower, color = condition), shape = 4, show.legend = F) + 
  geom_point(data = prediction_band, aes(x = time, y = upper, color = condition), shape = 4, show.legend = F) + 
  facet_wrap(~ name, scales = "free_y") +  # Facet by 'name' column
  labs(
    x = "Time",
    y = "Value",
    color = "Condition"
  ) +
  dMod::theme_dMod() +  # Apply dMod theme
  dMod::scale_color_dMod() +  # Apply dMod color scale to lines
  dMod::scale_fill_dMod()   # Apply the same color scale to the fill