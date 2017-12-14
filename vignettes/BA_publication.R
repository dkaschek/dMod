
library(dMod)
library(dplyr)
library(ggplot2)
theme_set(theme_dMod())

## Read the data ---------------------------------------------

data(BAdata)
data <- BAdata %>% 
  filter(experiment == "exp1") %>%
  as.datalist(split.by = c("experiment", "cations", "compound", "dose", "tca_time", "drug_time"))
covtable <- covariates(data)
print(covtable)

plot(data)


## Set up a base model ----------------------------------------

# Set up base reactions
reactions <- eqnlist() %>%
  addReaction("Tca_buffer", "Tca_cyto", "import_Tca*Tca_buffer", "Basolateral uptake") %>%
  addReaction("Tca_cyto", "Tca_buffer", "export_Tca_baso*Tca_cyto", "Basolateral efflux") %>%
  addReaction("Tca_cyto", "Tca_canalicular", "export_Tca_cana*Tca_cyto", "Canalicular efflux") %>%
  addReaction("Tca_canalicular", "Tca_buffer", "transport_Tca*Tca_canalicular", "Transport bile") %>%
  addReaction("0", "cations", "0", "Cation concentration in buffer")

# Translate to ODEs and let cation concentration change transport rate
odes <- reactions %>% 
  as.eqnvec() %>% 
  reparameterize("x ~ (x*crit/(crit + cations))", x = "transport_Tca")

# Set up events to control the experiment
events <- data.frame(
  var = c("Tca_buffer", "cations"),
  time = c(30, "time_remove_cations"),
  value = c(0, 0),
  method = c("replace", "replace")
)

# Make ODE model available as prediction function
x <- odemodel(odes, events = events, fixed = c("time_remove_cations", "cations"), modelname = "BA_basemodel") %>%
  Xs()

# Make observation function available
g <- eqnvec(
  buffer = "s*Tca_buffer",
  cellular = "s*(Tca_cyto + Tca_canalicular)") %>%
  Y(f = x, modelname = "obsfn", compile = TRUE)

# Parameterize the model
parameters <- getParameters(g, x)
transformation <- eqnvec() %>%
  reparameterize("x~x", x = parameters) %>%
  reparameterize("x~0", x = c("Tca_cyto", "Tca_canalicular")) %>%
  reparameterize("x~1", x = "cations")
  

p <- P()  
for (c in rownames(covtable)) {
  p <-p + transformation %>%
    reparameterize("x~cations", x = "time_remove_cations", cations = covtable[c, "cations"]) %>%
    reparameterize("x~exp(x)", x = parameters) %>%
    P(condition = c)
}

estimate <- getParameters(p)
pars <- structure(rep(-1, length(estimate)), names = estimate)
times <- seq(0, 200, 1)
(g*x*p)(times, pars) %>% plot(data = data, time <= 180)

obj <- normL2(data, g*x*p) + constraintL2(pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

(g*x*p)(times, myfit$argument) %>% plot(data = data)

profiles <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles %>% plotProfile(mode == "data")


# Add additional data
data$exp1_30_no_0_0_30 <- rbind(data$exp1_30_no_0_0_30, data.frame(
  name = "Tca_buffer", time = 0, value = 21.5, sigma = 2
))

obj <- normL2(data, g*x*p) + constraintL2(pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

(g*x*p)(times, myfit$argument) %>% plot(data = data)

profiles <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles %>% plotProfile(mode == "data")

# Explore parameter space
fits <- mstrust(obj, pars, rinit = 1, rmax = 10, samplefun = "runif", min = -5, max = 5, cores = 4)
parframe <- as.parframe(fits)
subframe <- unique(parframe[1:19,])

prediction <- predict(g*x*p, times = times, pars = subframe, data = data)

ggplot(prediction, aes(x = time, y = value, color = as.factor(.value))) + facet_wrap(~name*cations, scales = "free") +
  geom_line() +
  geom_point(data = as.data.frame(data), color = "black")

profiles <- lapply(1:nrow(subframe), function(i) {
  profile(obj, as.parvec(subframe, i), names(myfit$argument), cores = 4, algoControl = list(gamma = 10))
})
names(profiles) <- c("obj = 23", "obj = 24", "obj = 27")
profiles %>% plotProfile(mode == "data")

# Save objects from this session
save(g, x, p, transformation, data, file = "part1.RData")

  