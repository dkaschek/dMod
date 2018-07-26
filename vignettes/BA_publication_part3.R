
library(dMod)
library(dplyr)
library(ggplot2)
theme_set(theme_dMod())
rm(list = ls())

## Load contents of session 1 --------------------------------

load("part2.RData")

## Read the data ---------------------------------------------

data(BAdata)
data <- BAdata %>% 
  filter(experiment %in% c("exp2", "exp4") & compound %in% c("Csa", "Cpz")) %>%
  as.datalist()
covtable <- covariates(data)
print(covtable)

conditions_Csa <- rownames(covtable)[grepl("Csa", rownames(covtable))]
conditions_Cpz <- rownames(covtable)[grepl("Cpz", rownames(covtable))]

plot(data)

## Extend reactions by drug-transporter interaction

reactionsPD <- reactions %>%
  # Implement ability to return to original level 
  addReaction("0", "Bsep", "alpha_Bsep", "Transporter levels") %>%
  addReaction("Bsep", "0", "alpha_Bsep * Bsep", "Transporter levels") %>%
  addReaction("0", "Ntcp", "alpha_Ntcp", "Transporter levels") %>%
  addReaction("Ntcp", "0", "alpha_Ntcp * Ntcp", "Transporter levels") %>%
  addReaction("0", "Mrp3", "alpha_Mrp3", "Transporter levels") %>%
  addReaction("Mrp3", "0", "alpha_Mrp3 * Mrp3", "Transporter levels") %>%
  # Implement drug-transporter interaction
  addReaction("Bsep + Drug", "DrugBsep", "build_mech1 * Drug * Bsep", "DRUG effect") %>%
  addReaction("DrugBsep", "Bsep + Drug", "decay_mech1 * DrugBsep", "DRUG effect") %>%
  addReaction("Ntcp + Drug", "DrugNtcp", "build_mech1 * Drug * Ntcp", "DRUG effect") %>%
  addReaction("DrugNtcp", "Ntcp + Drug", "decay_mech1 * DrugNtcp", "DRUG effect") %>%
  addReaction("Mrp3 + Drug", "DrugMrp3", "build_mech2 * Drug * Mrp3", "DRUG effect") %>%
  addReaction("DrugMrp3", "Mrp3 + Drug", "decay_mech2 * DrugMrp3", "DRUG effect")

# Translate to ODEs and let cation concentration change transport rate
odesPD <- reactionsPD %>% 
  as.eqnvec() %>% 
  insert("transport_Tca ~ (transport_Tca * crit / (crit + cations))", rateconst = "transport_Tca") %>%
  insert("rate ~ (rate * transporter)", 
                 rate   = c("import_Tca", "export_Tca_baso", "export_Tca_cana"),
                 transporter = c("Ntcp" , "Mrp3"           , "Bsep"           ))

# Dosing events
eventsPD <- events %>%
  addEvent("Drug", "t_addDrug", "amount_Drug") %>%
  addEvent("Drug", "t_removeDrug", "0")

# Make ODE model available as prediction function
noSens <- parameters
xPD <- odemodel(odesPD, events = eventsPD, fixed = noSens, modelname = "BA_PDmodel") %>%
  Xs()

# Make observation function available
gPD <- eqnvec(
  buffer = "s*Tca_buffer",
  cellular = "s*(Tca_cyto + Tca_canalicular)") %>%
  Y(f = xPD, attach.input = FALSE, modelname = "obsfnPD", compile = TRUE)


# Parameterize the model
parameters <- getParameters(gPD, xPD)
pPD <- eqnvec() %>%
  define("x ~ x", x = parameters) %>%
  define("x ~ 0", x = c("Tca_cyto", "Tca_canalicular", "Tca_buffer", "Drug",
                      "DrugBsep", "DrugMrp3", "DrugNtcp")) %>%
  define("x ~ 1", x = c("cations", 
                      "Ntcp", "Bsep", "Mrp3")) %>%
  branch(covtable) %>% 
  define("t_removeCa ~ time", time = ifelse(changeCa == "yes", 30, 1e3)) %>%
  define("t_addTca ~ time", time = switch(as.character(experiment), exp2 = 0, exp4 = 30)) %>%
  define("t_removeTca ~ time", time = switch(as.character(experiment), exp2 = 30, exp4 = 1e3)) %>%
  define("t_addDrug ~ time", time = switch(as.character(experiment), exp2 = 30, exp4 = 0)) %>%
  define("t_removeDrug ~ time", time = switch(as.character(experiment), exp2 = 1e3, exp4 = 30)) %>%
  define("amount_Drug ~ dose*exp(DRUGSCALE)", dose = dose) %>%
  insert("x ~ exp(X)", x = parameters, X = toupper(parameters)) %>%
  # Make build parameters specific for compound
  insert("BUILD_x ~ BUILD_x_compound", 
                 x = c("MECH1", "MECH2"), 
                 compound = compound) %>%
  # Make decay parameter specific for compound
  insert("DECAY_x ~ DECAY_x_compound", 
                 x = c("MECH1", "MECH2"), 
                 compound = compound) %>%
  # Insert estimated values
  insert("value ~ estimate", 
         value = as.character(partable$name), 
         estimate = partable$value) %>%
  P()
  
  
estimate <- getParameters(pPD)

# pars <- log(c(
#   
#   ALPHA_BSEP = 1e-3,
#   ALPHA_MRP3 = 0.1,
#   ALPHA_NTCP = 1e-3, 
#   
#   BUILD_MECH1_Csa = 0.6,
#   BUILD_MECH1_Cpz = 0.001,
#   BUILD_MECH2_Csa = 0.02,
#   BUILD_MECH2_Cpz = 1e-9,
#   
#   DECAY_MECH1_Csa = 0.5,
#   DECAY_MECH1_Cpz = 1e-9,
#   DECAY_MECH2_Csa = 0.03,
#   DECAY_MECH2_Cpz = 1e-9,
#   
#   DRUGSCALE = 2.3
#   
# ))

times <- seq(0, 200, 1)
pars <- structure(rep(-10, length(estimate)), names = estimate)
(gPD*xPD*pPD)(times, pars, deriv = FALSE) %>% plot(data = data, time > 30 & grepl("exp2", condition))
(gPD*xPD*pPD)(times, pars, deriv = FALSE) %>% plot(data = data, time > 30 & grepl("exp4", condition))


obj <- normL2(data, gPD*xPD*pPD) + constraintL2(pars, sigma = 10)
pars <- structure(rep(-1, length(estimate)), names = estimate)
fits <- mstrust(obj, pars, fits = 50, rinit = 1, rmax = 10, 
                samplefun = "runif", min = -2, max = 2, cores = 4,
                conditions = conditions_Cpz) %>% as.parframe()
myfit <- trust(obj, pars, rinit = 1, rmax = 10)
(gPD*xPD*pPD)(times, myfit$argument) %>% plot(data = data, facet = "grid", time > 30)
(gPD*xPD*pPD)(times, as.parvec(fits)) %>% plot(data = data, facet = "grid")

plotValues(fits, .1, value < 1000)
plotPars(fits, .1, value < 1000)

subframe <- fits[fits$value < 1e3 & fits$converged, ] %>% unique()

pred <- predict(gPD*xPD*pPD, times = times, pars = subframe, data = data)
ggplot(pred, aes(x = time, y = value, color = as.factor(.value))) + 
  facet_wrap( ~ condition*name, scales = "free") +
  geom_point(data = attr(pred, "data"), color = "black") + geom_line()

pred <- predict(xPD*pPD, times = times, pars = subframe, data = data)
ggplot(pred, aes(x = time, y = value, color = as.factor(.value))) + 
  facet_wrap( ~ condition*name, scales = "free") + geom_line()

bestfit <- as.parvec(fits)
profiles_part3 <- profile(obj, bestfit, names(bestfit), cores = 4, limits = c(-10, 10), 
                          stepControl = list(stop = "data"),
                          conditions = conditions_Cpz)
profiles_part3 %>% plotProfile(mode == "data")

mypars <- bestfit[setdiff(names(bestfit), "DRUGSCALE"), drop = TRUE]
myfixed <- bestfit["DRUGSCALE", drop = TRUE]
profiles_part3 <- profile(obj, mypars, names(mypars), cores = 4, limits = c(-10, 10), 
                          stepControl = list(stop = "data"),
                          conditions = conditions_Cpz,
                          fixed = myfixed)
profiles_part3 %>% plotProfile(mode == "data")
