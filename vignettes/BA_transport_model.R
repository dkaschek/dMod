##
## Replication script to reproduce the results shown in the manuscript
##


library(dMod)
library(ggplot2)
theme_set(theme_dMod())
setwd(tempdir())

# Section 4.1. Simulation and prediction ------------------------------

# Add reactions
reactions <- NULL
reactions <- addReaction(reactions, "TCA_buffer", "TCA_cell",
                         rate = "import*TCA_buffer",
                         description = "Uptake")
reactions <- addReaction(reactions, "TCA_cell", "TCA_buffer",
                         rate = "export_sinus*TCA_cell",
                         description = "Sinusoidal export")
reactions <- addReaction(reactions, "TCA_cell", "TCA_cana",
                         rate = "export_cana*TCA_cell",
                         description = "Canalicular export")
reactions <- addReaction(reactions, "TCA_cana", "TCA_buffer",
                         rate = "reflux*TCA_cana",
                         description = "Reflux into the buffer")
# Translate into ODE model
mymodel <- odemodel(reactions, modelname = "bamodel")

# Generate prediction function from ODE model
x <- Xs(mymodel, condition = NULL)

times <- seq(0, 50, .1)
pars <- c(TCA_buffer = 1,
          TCA_cell = 0,
          TCA_cana = 0,
          import = 0.2,
          export_sinus = 0.2,
          export_cana = 0.04,
          reflux = 0.1)

# Generate the prediction
out <- x(times, pars)
plot(out)
# Get sensitivities for the rate parameters only, pars[4:7]
out <- getDerivs(x(times, pars[4:7], fixed = pars[1:3]))
plot(out)


# Section 4.2. Observation function and simulated data ---------------

# Generate observation function
observables <- eqnvec(
  buffer = "s*TCA_buffer",
  cellular = "s*(TCA_cana + TCA_cell)"
)
g <- Y(observables, x, condition = NULL,
       compile = TRUE, modelname = "obsfn")

# Reset parameter values
pars["TCA_cell"] <- 0.3846154
pars["TCA_cana"] <- 0.1538462
pars["TCA_buffer"] <- 0
pars["s"] <- 1e3

out <- (g*x)(times, pars, conditions = "standard")

# Simuate data
set.seed(1)
timesD <- c(0.1, 1, 3, 7, 11, 15, 20, 41)
datasheet <- subset(as.data.frame(out),
                    time %in% timesD & name %in% names(observables))

datasheet$sigma <- sqrt(datasheet$value + 1)
datasheet$value <- rnorm(nrow(datasheet), datasheet$value, datasheet$sigma)
data <- as.datalist(datasheet)
plot(out, data)

# Section 4.3. Parameter transformation ------------------------------

p <- P(
  trafo = eqnvec(
    TCA_buffer = "0",
    TCA_cell = "exp(TCA_cell)",
    TCA_cana = "exp(TCA_cana)",
    import = "exp(import)",
    export_sinus = "exp(export_sinus)",
    export_cana = "exp(export_cana)",
    reflux = "exp(reflux)",
    s = "exp(s)"
  ),
  condition = "standard"
)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*p)(times, pouter), data)

# 4.4. Objective function and model fitting --------------------------

obj <- normL2(data, g*x*p) + constraintL2(pouter, sigma = 10)
myfit <- trust(obj, pouter, rinit = 1, rmax = 10)
plot((g*x*p)(times, myfit$argument), data)

out_mstrust <- mstrust(obj, pouter, rinit = 1, rmax = 10, iterlim = 500,
                       sd = 4,
                       cores = 4, fits = 50)

myframe <- as.parframe(out_mstrust)
P <- plotValues(myframe, tol = .01, value < 100)
print(P)
plotPars(myframe, tol = .01, value < 100)

# Make predictions for the different log-likelihood values
select <- attr(P, "jumps")
prediction <- predict(g*x*p, times = times, pars = myframe[select,], data = data)
ggplot(prediction, aes(x = time, y = value, group = .index, color = .value)) + 
  facet_wrap(~name, scales = "free") + geom_line() +
  geom_point(data = attr(prediction, "data")) +
  theme_dMod()

# 4.5. Working with several conditions -------------------------------

pars["reflux"] <- 1e3
out <- (g*x)(times, pars, conditions = "open")
datasheet <- subset(as.data.frame(out),
                    time %in% timesD & name %in% names(observables))
datasheet$sigma <- sqrt(datasheet$value + 1)
datasheet$value <- rnorm(nrow(datasheet), datasheet$value, datasheet$sigma)
data <- data + as.datalist(datasheet)

trafo <- summary(p)$standard$equations
trafo["reflux"] <- "exp(reflux_open)"
p <- p + P(trafo, condition = "open")

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
obj <- normL2(data, g*x*p) + constraintL2(pouter, sigma = 10)

out_mstrust <- mstrust(obj, pouter, rinit = 1, rmax = 10, iterlim = 500,
                       sd = 4, cores = 4, fits = 50)
myframe <- as.parframe(out_mstrust)
plotValues(myframe, tol = 1, value < 100)
plotPars(myframe, tol = 1, value < 100)
bestfit <- as.parvec(myframe)
plot((g*x*p)(times, bestfit), data)

# 4.6. Parameter uncertainty and identifiability ---------------------

profiles <- profile(obj, bestfit, names(bestfit), limits = c(-5, 5), cores = 4, verbose = TRUE)
plotProfile(profiles, mode == "data")
plotPaths(profiles, whichPar = "s")

# 4.7. Steady-state constraints and implicit transformations ---------

pSS <- NULL
trafos <- summary(p)
conditions <- names(trafos)
for (n in conditions) {
  equations <- trafos[[n]]$equations
  equations["TCA_cana"] <- "exp(export_cana)*exp(TCA_cell)/exp(reflux)"
  pSS <- pSS + P(equations, condition = n)
}
outerpars <- getParameters(pSS)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

obj <- normL2(data, g*x*pSS) + constraintL2(pouter, sigma = 10)
bestfit <- trust(obj, pouter, rinit = 1, rmax = 10)$argument

profiles_SS_analytic <- profile(obj, bestfit, names(bestfit), limits = c(-5, 5), cores = 4)

plotProfile(list(noSS = profiles, SS_explicit = profiles_SS_analytic), mode == "data")


# Add reactions
reactions <- NULL
reactions <- addReaction(reactions, "TCA_buffer", "TCA_cell",
                         rate = "import*TCA_buffer",
                         description = "Uptake")
reactions <- addReaction(reactions, "TCA_cell", "TCA_buffer",
                         rate = "export_sinus*TCA_cell",
                         description = "Sinusoidal export")
reactions <- addReaction(reactions, "TCA_cell", "TCA_cana",
                         rate = "export_cana*TCA_cell",
                         description = "Canalicular export")
reactions <- addReaction(reactions, "TCA_cana", "TCA_buffer",
                         rate = "(reflux*(1-switch) + reflux_open*switch)*TCA_cana",
                         description = "Reflux into the buffer")
reactions <- addReaction(reactions, "0", "switch",
                         rate = "0",
                         description = "Create a switch")
# Translate into ODE model
mymodel <- odemodel(reactions, modelname = "bamodel")

# Set up implicit parameter transformation
f <- as.eqnvec(reactions)[c("TCA_buffer", "TCA_cana", "TCA_cell")]

f["TCA_cell"] <- "TCA_buffer + TCA_cana + TCA_cell - TCA_tot"
pSS <- P(f, method = "implicit",
         compile = TRUE, modelname = "pfn")

# Set up explicit parameter transformation
innerpars <- unique(c(getParameters(mymodel),
                      getSymbols(observables),
                      getSymbols(f)))

trafo <- repar("x ~ x", x = innerpars)
trafo <- repar("x ~ 0", x = reactions$states, trafo)
trafo <- repar("x ~ exp(x)", x = innerpars, trafo)

p <- P(trafo)

# Set up prediction function with events
event.buffer <- data.frame(var = "TCA_buffer",
                           time = 0,
                           value = 0,
                           method = "replace")
event.open <- data.frame(var = "switch",
                         time = 0,
                         value = 1,
                         method = "replace")
x <- Xs(mymodel,
        events = event.buffer,
        condition = "standard") +
  Xs(mymodel,
     events = rbind(event.buffer, event.open),
     condition = "open")

# Generate observation function with modified states/parameters
g <- Y(observables, x,
       compile = TRUE, modelname = "obsfn")

# Generate objective function
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
obj <- normL2(data, g*x*pSS*p) + constraintL2(pouter, sigma = 10)

bestfit <- trust(obj, pouter, rinit = 1, rmax = 10)$argument

profiles_SS_implicit <- profile(obj, bestfit, names(bestfit), limits = c(-5, 5), cores = 4)

plotProfile(list(noSS = profiles, 
                 SS_explicit = profiles_SS_analytic, 
                 SS_implicit = profiles_SS_implicit),
            mode == "data")


# 4.8. Prediction uncertainty and validation profiles ----------------

obj.validation <- datapointL2(name = "TCA_cell",
                              time = 41,
                              value = "d1",
                              sigma = .1,
                              condition = "standard")

myfit <- trust(obj + obj.validation,
               parinit = c(d1 = 1, bestfit[-7]),
               fixed = c(TCA_tot = bestfit[[7]]),
               rinit = 1, rmax = 10)

profile_prediction <- profile(obj + obj.validation,
                              myfit$argument, "d1", limits = c(-5, 5),
                              stepControl = list(stop = "data"),
                              fixed = c(TCA_tot = bestfit[[7]]))

plotProfile(profile_prediction, mode %in% c("data", "validation"))


# 4.8.1 Prediction band (prediction uncertainty for several time points) --------------
prediction_band <- do.call(rbind, lapply(c(7., 11., 20., 41.), function(t) {
  
  cat("Computing prediction profile for t =", t, "\n")
  
  obj.validation <- datapointL2(name = "TCA_cell",
                                time = t,
                                value = "d1",
                                sigma = .1,
                                condition = "standard")
  
  myfit <- trust(obj + obj.validation,
                 parinit = c(d1 = 1, bestfit[-7]),
                 fixed = c(TCA_tot = bestfit[[7]]),
                 rinit = 1, rmax = 10)
  
  profile_prediction <- profile(obj + obj.validation,
                                myfit$argument, "d1", limits = c(-5, 5),
                                fixed = c(TCA_tot = bestfit[[7]]))
  
  d1 <- confint(profile_prediction, val.column = "value")
  
  # Output
  data.frame(time = t, condition = "standard", name = "TCA_cell",  d1[-1])
  
  
}))

plot((g*x*pSS*p)(times, bestfit), data) + 
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), 
              data = prediction_band,
              lty = 0, alpha = .3)
