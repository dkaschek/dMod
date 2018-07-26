

# Parameter Estimation in an ODE Model of Bile Acid Transport

## Preliminaries

library(dMod)
library(dplyr)
library(ggplot2)
library(pander)
theme_set(theme_dMod())

setwd("/tmp")


## The model scheme

## Read the data


data(BAdata)
data <- BAdata %>% 
  filter(experiment == "exp1") %>%
  as.datalist()
covtable <- covariates(data)
plot(data)

## Setting up the base model

### Setting up the prediction function


reactions <- eqnlist() %>%
 #addReaction(from, to, rate, description) %>%  
  addReaction("Tca_buffer", "Tca_cyto", "import_Tca*Tca_buffer", "Basolateral uptake") %>%
  addReaction("Tca_cyto", "Tca_buffer", "export_Tca_baso*Tca_cyto", "Basolateral efflux") %>%
  addReaction("Tca_cyto", "Tca_canalicular", "export_Tca_cana*Tca_cyto", "Canalicular efflux") %>%
  addReaction("Tca_canalicular", "Tca_buffer", "transport_Tca*Tca_canalicular", "Transport bile") %>%
  addReaction("0", "cations", "0", "Cation concentration in buffer") 

odes <- reactions %>% 
  as.eqnvec() %>% 
  insert("transport_Tca ~ (transport_Tca*crit/(crit + cations))")

events <- eventlist() %>%
  #addEvent(var, time, value, method = "replace") %>%
  addEvent("Tca_buffer", "t_addTca", "amount_Tca") %>%
  addEvent("Tca_buffer", "t_removeTca", "0") %>%
  addEvent("cations", "t_removeCa", "0")

noSens <- c("t_removeCa", "t_addTca", "t_removeTca", "cations")
x <- odemodel(odes, events = events, fixed = noSens, modelname = "BA_basemodel") %>%
  Xs()

## Demonstration of the usage of x()
n <- getParameters(x)
pars <- runif(length(n), min = 0, max = 100)
names(pars) <- n

times <- 0:200
pred <- x(times, pars)
plot(pred) + theme(legend.position = "none")

### Setting up the observation function

g <- eqnvec(buffer = "s*Tca_buffer",
            cellular = "s*(Tca_cyto + Tca_canalicular)") %>%
  Y(f = x, modelname = "obsfn", compile = TRUE)

# Demonstration of the usage of g
pars["s"] <- 1
pred <- (g*x)(times, pars)
plot(pred) + theme(legend.position = "none")


### Parameterization



parameters <- getParameters(g, x)
p <- eqnvec() %>%
  define("x ~ x", x = parameters) %>%
  define("x ~ 0", x = c("Tca_cyto", "Tca_canalicular", "Tca_buffer")) %>%
  define("cations ~ 1") %>%
  define("t_addTca ~ 0") %>%
  define("t_removeTca ~ 30") %>%
  branch(covtable) %>% 
  define("t_removeCa ~ time", time = ifelse(changeCa == "yes", 30, 1e3)) %>%
  insert("x ~ exp(X)", x = parameters, X = toupper(parameters)) %>%
  P()


## Parameter estimation

set.seed(12837)
estimate <- getParameters(p)
pars <- structure(runif(length(estimate), min = -2, max = 0), names = estimate)
times <- seq(0, 200, 1)
(g*x*p)(times, pars) %>% plot(data = data, time <= 180)


obj <- normL2(data, g*x*p) + constraintL2(0*pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

out <- mstrust(obj, pars, rinit = 1, rmax = 10, cores = 4) %>% as.parframe()

(g*x*p)(times, myfit$argument) %>% plot(data = data, time <= 180)
### Identifiability
profiles <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles %>% plotProfile(mode == "data")


# Add additional data
data$exp1_yes_no_0 <- data$exp1_yes_no_0 %>%
  rbind(data.frame(name = "Tca_buffer", time = 1, value = 21.5, sigma = 2))

obj <- normL2(data, g*x*p) + constraintL2(0*pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

pred <- (g*x*p)(times, myfit$argument) 
pred %>% plot(data = data, name %in% c("buffer", "cellular", "Tca_buffer"))


profiles <- profile(obj, myfit$argument, c("S", "AMOUNT_TCA"), cores = 2)
profiles %>% plotProfile(mode == "data")


# Explore parameter space
fits <- mstrust(obj, pars, fits = 30, rinit = 1, rmax = 10, samplefun = "runif", min = -5, max = 5, cores = 4)
parframe <- as.parframe(fits)
plotValues(parframe[1:28,], , value < 100)
plotPars(parframe[1:28,], , value < 100)
subframe <- unique(parframe[1:20,])

prediction <- predict(g*x*p, times = times, pars = subframe, data = data)

ggplot(prediction, aes(x = time, y = value, color = as.factor(round(.value)))) + facet_wrap( ~ name*changeCa, scales = "free") +
  geom_line() +
  geom_point(data = as.data.frame(data), color = "black") +
  scale_color_discrete(name = "obj. value")


profiles <- lapply(1:nrow(subframe), function(i) {
  profile(obj, as.parvec(subframe, i), names(myfit$argument), limits = c(-4, 4), cores = 4)
})
names(profiles) <- c("obj = 23", "obj = 24", "obj = 27")
profiles %>% plotProfile(mode == "data")


partable <- confint(profiles[[1]], level = 0.9, val.column = "value")


# Save objects from this session
save.image(file = "part1.RData")

  

