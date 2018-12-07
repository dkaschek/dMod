# publised in https://www.ncbi.nlm.nih.gov/pubmed/27811075
# D2D: transferred to D2D exmple models, named bruno
# simple model, sort of A -> B -> C
# no error model fitting

## Model Definition ------------------------------------------------------

# Read in model csv
reactionlist <- read.csv(paste0(system.file(package = "dMod"),'/examples/example_CCD4/model.csv')) 

f <- as.eqnlist(reactionlist)

# Define new observables based on ODE states
species <- f$states

observables <- c(
  obcar = "bcar",
  obcry = "bcry",
  obio = "bio",
  ob10 = "b10",
  oohb10 = "ohb10",
  ozea = "zea"
)

# List of fixed parameters which are known beforehand
fixed <- NULL
forcings <- NULL 

# Add observable ODEs to the original ODEs or use an observation function
g <- Y(observables, as.eqnvec(f), compile = TRUE,modelname = "obs", attach.input = FALSE)

# Generate the model C files, compile them and return a list with func and extended.
do.compile <- TRUE
cat("now compile the model", do.compile,"\n")
model0 <- odemodel(as.eqnvec(f), fixed = fixed, forcings = NULL, jacobian = "inz.lsodes", compile = do.compile, modelname = "ccd4")

## Data ----------------------------------------------------------------------

# Data was preprocessed using an error model fit
# datasheet <- subset(read.table("data_mark_out.csv",sep=",",header=TRUE), time >0 ) # needs columns condition, name, time, value, sigma, time = 0 were removed due to experimental problems
# dMod::fitErrorModel(datasheet, factors = c("name"), plotting = FALSE, errorModel = "(exp(s0) + x^2*exp(srel))", par = c(s0 = -6, srel = -4))$sigma

# Final data used for the publication (class: datalist):
data <- readRDS(paste0(system.file(package = "dMod"),'/examples/example_CCD4/datalist.Rds'))
data <- as.datalist(data)
conditions <- getConditions(data)
## Parameter Transformations -------------------------------------------

# Define inner parameters (parameters occurring in the equations except forcings)
# Add names(observables) if addObservables(observables, f) is used
innerpars <- getSymbols(c(as.character(as.eqnvec(f)), names(as.eqnvec(f)), as.character(observables)), exclude=c(forcings, "time"))
names(innerpars) <- innerpars

# Define additional parameter constraints, e.g. steady-state conditions
# Parameters (left-hand side) are replaced in the right-hand side of consecutive lines by resolveRecurrence() 
# setting the initial value to zero for most states
constraints <- c(
  b10 = 0,
  bio = 0,
  ohb10 = 0,
  ohbio = 0,
  zea = 0
)
constraints <- resolveRecurrence(constraints)

# Build up a parameter transformation (constraints, log-transform, etc.)
# Start with replacing initial value parameters of the observables
trafo <- replaceSymbols(names(observables), observables, innerpars)
# Then employ the other parameter constraints
trafo <- replaceSymbols(names(constraints), constraints, trafo)
# Then do a log-transform of all parameters (if defined as positive numbers)
parsT <- unique(c(getSymbols(constraints), innerpars))
trafo <- replaceSymbols(parsT, paste0("exp(log", parsT, ")"), trafo)

parameters <- innerpars[!(innerpars %in% species)]

offsets <- paste0("log",parameters[grepl("off", parameters) & !(grepl("zea", parameters))& !(grepl("bcry", parameters))])

## Specify different conditions -----------------------------------------------------
# Set condition-specific parameter transformations
# Mostly the conditions differ by the set of initial values for the states (different substrates provided for each condition)
trafoLTot <- lapply(conditions, function(con){ 
  trafo
  if(grepl("zea", con)){
    # the conditions betacar_zea and zea were not measured in parallel with the other conditions
    specific_zea <- c(parameters[grepl("k", parameters)])
    trafo[specific_zea] <- paste0(trafo[specific_zea], "*exp(logk_zea)")
  }
  if(grepl("betacry",con)){
     trafo["kb1"] <- "0"
     trafo["bcar"] <- "0"
  }
  if(grepl("betacar", con)){
    trafo["bcry"] <- "0"
    trafo["kc1"] <- "0"
    trafo["kc2"] <- "0"
    trafo["kc4"] <- "0"
  }
  if(con == "betacar_zea"){
    trafo["bcar"] <- "exp(logbcar_zea)"
  }
  if(con == "zea"){
    trafo["zea"] <- "exp(logzeaini)"
    trafo["bcry"] <- "0"
    trafo["bcar"] <- "0"
    trafo["kc1"] <- "0"
    trafo["kc2"] <- "0"
    trafo["kb1"] <- "0"
    trafo["kb2"] <- "0"
    trafo["off_ob10"] <- "0"
  }
  if(con == "b10"){
    trafo["b10"] <- "exp(logb10ini)"
    trafo["bcar"] <- "0"
    trafo["bcry"] <- "0"
    trafo["kb1"] <- "0"
    trafo["kc1"] <- "0"
    trafo["kc2"] <- "0"
    trafo["k5"] <- "0"
    trafo["off_oohb10"] <- "0"
  }  
  if(con == "ohb10"){
    trafo["ohb10"] <- "exp(logohb10ini)"
    trafo["bcar"] <- "0"
    trafo["bcry"] <- "0"
    trafo["kb1"] <- "0"
    trafo["kc1"] <- "0"
    trafo["kc2"] <- "0"
    trafo["k5"] <- "0"
    trafo["off_ob10"] <- "0"
  }  
  if(con == "betacar_b10"){
    trafo["b10"] <- "exp(logbcar)"
    trafo["bcar"] <- "0"
    trafo["bcry"] <- "0"
    trafo["kc1"] <- "0"
    trafo["kc2"] <- "0"
    trafo["kc4"] <- "0"
  }
    return(trafo)
}) 
names(trafoLTot) <- conditions

p <- NULL
for (C in conditions) {
  p <- p + P(trafoLTot[[C]], condition = C)
}

optionsOde = list(method = "lsoda")#, lrw = 524)
optionsSens = list(method = "lsodes", rtol = 1e-8, atol = 1e-8)
x <- Xs(model0, optionsOde = optionsOde, optionsSens = optionsSens)

trafoL <- trafoLTot[conditions]

## Objective Function -------------------------------------------------------

# Initalize parameters 
outerpars <- getSymbols(do.call(c, trafoL[conditions]))
set.seed(2)
# building a weak parameter prior that can be used to ease optimization
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
pouter <- rnorm(length(prior), prior, 1); names(pouter) <- outerpars
cOuter <- constraintL2(mu = prior[names(pouter)], sigma = 10)

# Objective function for trust() 
obj <- normL2(data, g*x*p) # + cOuter

