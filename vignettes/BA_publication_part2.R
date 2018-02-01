
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

parameters <- getParameters(g, x)
p_add <- getEquations(p, conditions = "exp1_no_no_0") %>%
  define("t_addTca ~ 60") %>%
  define("x ~ 1e3", x = c("t_removeTca", "t_removeCa")) %>%
  P(condition = "exp3_no_no_0")

estimate <- getParameters(p, p_add)
pars <- structure(rep(-1, length(estimate)), names = estimate)
pars[partable$name] <- partable$value

times <- seq(0, 200, 1)
(g*x*(p + p_add))(times, pars) %>% plot(data = data, time <= 180)

y <- g * x * (p + p_add)
obj <- normL2(data, y) + constraintL2(pars, sigma = 10)
myfit <- trust(obj, pars, rinit = 1, rmax = 10, iterlim = 500)

y(times, myfit$argument) %>% plot(data = data)

profiles_part2 <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles_part2 %>% plotProfile(mode == "data")

plotProfile(list(old = profiles[[1]], new = profiles_part2), mode == "data")

partable <- confint(profiles_part2, level = 0.68)

# Save objects from this session
save.image(file = "part2.RData")

  