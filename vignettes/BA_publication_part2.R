
library(dMod)
library(dplyr)
library(ggplot2)
theme_set(theme_dMod())
rm(list = ls())

## Load contents of session 1 --------------------------------

load("part1.RData")
loadDLL(g, x)

## Read the data ---------------------------------------------

data(BAdata)
data <- data + BAdata %>% 
  filter(experiment == "exp3" & compound == "no") %>%
  as.datalist()
covtable <- covariates(data)
print(covtable)

plot(data)


## Add new condition to parameterization --------------------------


# Parameterize the model
parameters <- getParameters(g, x)
transformation <- eqnvec() %>%
  reparameterize("x~x", x = parameters) %>%
  reparameterize("x~0", x = c("Tca_cyto", "Tca_canalicular", "Tca_buffer", "value_cations")) %>%
  reparameterize("x~1", x = "cations")
  

for (c in rownames(covtable)[3]) {
  p <-p + transformation %>%
    reparameterize("x~cations", x = "change_cations", cations = covtable[c, "cations"]) %>%
    reparameterize("x~buffer", x = "change_buffer", buffer = covtable[c, "tca_time"]) %>%
    reparameterize("x~exp(x)", x = parameters) %>%
    P(condition = c)
}

estimate <- getParameters(p)
pars <- structure(rep(-1, length(estimate)), names = estimate)
pars[partable$name] <- partable$value

times <- seq(0, 200, 1)
(g*x*p)(times, pars) %>% plot(data = data, time <= 180)

obj <- normL2(data, g*x*p) + constraintL2(pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

(g*x*p)(times, myfit$argument) %>% plot(data = data)

profiles_part2 <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles_part2 %>% plotProfile(mode == "data")

plotProfile(list(old = profiles[[1]], new = profiles_part2), mode == "data")

partable <- confint(profiles_part2, level = 0.68)

# Save objects from this session
save.image(file = "part2.RData")

  