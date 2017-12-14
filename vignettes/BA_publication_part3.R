
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
  filter(experiment == "exp3" & compound != "no") %>%
  as.datalist()
covtable <- covariates(data)
print(covtable)

plot(data)

## Extend reactions by drug-transporter interaction

reactions_withDrug <- reactions %>%
  # Implement ability to return to original level 
  addReaction("0", "Bsep", "alpha_Bsep", "Transporter levels") %>%
  addReaction("Bsep", "0", "alpha_Bsep * Bsep", "Transporter levels") %>%
  addReaction("0", "Ntcp", "alpha_Ntcp", "Transporter levels") %>%
  addReaction("Ntcp", "0", "alpha_Ntcp * Ntcp", "Transporter levels") %>%
  addReaction("0", "Mrp3", "alpha_Mrp3", "Transporter levels") %>%
  addReaction("Mrp3", "0", "alpha_Mrp3 * Mrp3", "Transporter levels") %>%
  # Implement drug-transporter interaction
  addReaction("Bsep + Drug", "DrugBsep", "build_DrugBsep * Drug * Bsep", "DRUG effect") %>%
  addReaction("DrugBsep", "Bsep + Drug", "decay_DrugBsep * DrugBsep", "DRUG effect") %>%
  addReaction("Ntcp + Drug", "DrugNtcp", "build_DrugNtcp * Drug * Ntcp", "DRUG effect") %>%
  addReaction("DrugNtcp", "Ntcp + Drug", "decay_DrugNtcp * DrugNtcp", "DRUG effect") %>%
  addReaction("Mrp3 + Drug", "DrugMrp3", "build_DrugMrp3 * Drug * Mrp3", "DRUG effect") %>%
  addReaction("DrugMrp3", "Mrp3 + Drug", "decay_DrugMrp3 * DrugMrp3", "DRUG effect")

# Translate to ODEs and let cation concentration change transport rate
odes <- reactions_withDrug %>% 
  as.eqnvec() %>% 
  reparameterize("rateconst ~ (rateconst * crit / (crit + cations))", rateconst = "transport_Tca") %>%
  reparameterize("rateconst ~ (rateconst * transporter)", 
                 rateconst   = c("import_Tca", "export_Tca_baso", "export_Tca_cana"),
                 transporter = c("Ntcp"      , "Mrp3"           , "Bsep"           ))

 # Set up events to control the experiment
events <- data.frame(
  var = c("Tca_buffer", "cations", "Drug"),
  time = c("change_buffer", "change_cations", "change_drug"),
  value = c("value_buffer", "value_cations", "value_drug"),
  method = c("replace", "replace", "replace")
) 

# Make ODE model available as prediction function
noSens <- c("change_buffer", "change_cations", "value_cations", "cations", 
            "change_drug",
            "import_Tca", "export_Tca_baso", "export_Tca_cana", "transport_Tca", "crit", "value_buffer")
xPD <- odemodel(odes, events = events, fixed = noSens, modelname = "BA_PDmodel") %>%
  Xs()

# Make observation function available
gPD <- eqnvec(
  buffer = "s*Tca_buffer",
  cellular = "s*(Tca_cyto + Tca_canalicular)") %>%
  Y(f = xPD, modelname = "obsfnPD", compile = TRUE)


# Parameterize the model
parameters <- getParameters(g, xPD)
transformation <- eqnvec() %>%
  reparameterize("x~x", x = parameters) %>%
  reparameterize("x~0", x = c("Tca_cyto", "Tca_canalicular", "Tca_buffer", "value_cations", "DrugBsep", "DrugMrp3", "DrugNtcp")) %>%
  reparameterize("x~1", x = c("cations", "Ntcp", "Bsep", "Mrp3"))
  
p <- P()
for (c in rownames(covtable)) {
  p <-p + transformation %>%
    reparameterize("x~cations", x = "change_cations", cations = covtable[c, "cations"]) %>%
    reparameterize("x~buffer", x = c("change_buffer", "change_drug"), buffer = covtable[c, "tca_time"]) %>%
    reparameterize("x~amount", x = "Drug", amount = covtable[c, "dose"]) %>%
    reparameterize("x~drug", x = "value_drug", drug = "0") %>%
    reparameterize("x~exp(x)", x = parameters) %>%
    reparameterize("build_x~build_x_compound", x = c("DrugBsep", "DrugMrp3", "DrugNtcp"), compound = covtable[c, "compound"]) %>%
    reparameterize("decay_x~decay_x_compound", x = c("DrugBsep", "DrugMrp3", "DrugNtcp"), compound = covtable[c, "compound"]) %>%
    reparameterize("value~estimate", value = as.character(partable$name), estimate = partable$value) %>%
    P(condition = c)
}


estimate <- getParameters(p)
pars <- structure(rep(0, length(estimate)), names = estimate)

times <- seq(0, 200, 1)
(gPD*xPD*p)(times, pars) %>% plot(data = data, time <= 180)

obj <- normL2(data, gPD*xPD*p) + constraintL2(pars, sigma = 10)
fixed <- pars["drugscale"]
myfit <- trust(obj, pars[setdiff(names(pars), names(fixed))], rinit = 1, rmax = 10, iterlim = 100, fixed = fixed)

(gPD*xPD*p)(times, myfit$argument, fixed = fixed) %>% plot(data = data)

profiles_part2 <- profile(obj, myfit$argument, names(myfit$argument), cores = 4)
profiles_part2 %>% plotProfile(mode == "data")


