conditions <-c(
  #"pfl_parB", # mutant B: would need specific parameters
  "standard_par",
  #"standard_parB", # mutant B: would need specific parameters
  "standard_parh",
  "pfl_par"
)

## Model Definition ------------------------------------------------------

nmodel <- "simple_pfcoop" #pub model
#nmodel <- "simple_pfc_nopfc" # pub ohne coop

# Read in model csv
reactionlist <- read.csv(paste0("inst/examples/example_phytoene/phyto_",nmodel,".csv") )
if(nmodel == "simple_pfcoop")
  reactionlist <- reactionlist[c(1,3,5,7,8),]

# Translate data.frame into equations
f <- as.eqnlist(reactionlist)

# add equations for enzyme degradation
f <- addReaction(f, "Fad", "Fadage", rate="amount_D*Fad")
f <- addReaction(f, "Fadr", "Fadage", rate="amount_D*Fadr")

species <- f$states

errors <- NULL

observables <- c(
    op = "s1*(p)",
    opf = "s2*(pf+pfc)",
    oz = "s3*z",
    ptot_obs ="s1*(p)+s2*(pf+pfc)+s3*z", 
    ptot ="p+pf+z+pfc",
    fadtot="Fad+Fadr"
)

if(do.fiterrors){
  errors <- as.eqnvec(
    c(
      paste0("op*sigma_op_rel"),
      paste0("opf*sigma_opf_rel + sigma_opf_abs"),
      paste0("oz*sigma_oz_rel + sigma_oz_abs")
      ),
    names = names(observables)[1:3]
  )[1:3]
  }
  

# List of fixed parameters which are known beforehand
fixed <- NULL
forcings <- NULL 

# Add observable ODEs to the original ODEs or use an observation function
# Choose one of the three options, or combine them
g <- Y(observables, as.eqnvec(f), modelname = "obs_mono", compile=TRUE)
err <- NULL
if (!is.null(errors)) err <- Y(errors, c(observables, as.eqnvec(f)), states = names(errors), 
                               attach.input = FALSE, modelname = "err", compile = TRUE)


# Generate the model C files, compile them and return a list with func and extended.
cat("now compile the model", do.compile,"\n")
model0 <- odemodel(as.eqnvec(f), fixed = fixed, forcings = NULL, jacobian = "inz.lsodes", compile = do.compile,modelname = "mono")

## Parameter Transformations -------------------------------------------

constraints <- c(
  Fadage = "0", # dead enzyme
  Fad = "0.18", # given/measured by biologists
  Fadr="0",
  sigma_opf_abs = "0",
  DPQ="17.5*0.55",
  s1 = "1",
  s2 = "(l1/l2)*s1",
  s3 = "l2*s2",
  # standard setup/condition:
  pf = "0", 
  pfc = "0",
  z="0"
)

constraints <- resolveRecurrence(constraints)

# Define inner parameters (parameters occurring in the equations except forcings)
# Add names(observables) if addObservables(observables, f) is used
innerpars <- getSymbols(c(as.character(as.eqnvec(f)), names(as.eqnvec(f)), as.character(observables), constraints, errors), exclude=c(forcings, "time", names(observables)))
names(innerpars) <- innerpars

# Build up a parameter transformation (constraints, log-transform, etc.)
# Start with replacing initial value parameters of the observables
trafo <- replaceSymbols(names(observables), observables, innerpars)
# Then employ the other parameter constraints
trafo <- replaceSymbols(names(constraints), constraints, trafo)
# Then do a log-transform of all parameters (if defined as positive numbers)
logpars <- innerpars
trafo <- replaceSymbols(logpars, paste0("exp(log", logpars, ")"), trafo)

## Specify different conditions -----------------------------------------------------

# Set condition-specific parameter transformations 
specific <- c("logp")

trafoLTot <- lapply(conditions, function(con){ 
  trafo
  if(con %in% c("standard_par","standard_parB","standard_parh","pfl_par","pfl_parB")){
    specific_new <- c("logl1", "logl2")
    trafo <- replaceSymbols(specific_new, paste(specific_new, "new", sep="_"), trafo)
  }
  if(con %in% c("standard_par","standard_parB","standard_parh")){
    trafo <- replaceSymbols(specific, paste(specific, con, sep="_"), trafo)
    }
  if(con %in% c("pfl_par","pfl_parB")){
    trafo["p"] <- 0
    trafo <- replaceSymbols("logp", "0", trafo)
    trafo["kooperpFad"] <- "0"
    trafo["pf"] <- paste0("exp(logpf_", con,")")
  }
  return(trafo)
}) 
names(trafoLTot) <- conditions

trafoL <- trafoLTot[conditions]
p <- x <- NULL
optionsOde = list(method = "lsoda")
optionsSens = list(method = "lsodes", rtol = 1e-8, atol = 1e-8)

for (C in conditions) {
  p <- p + P(trafoLTot[[C]], condition = C)
  events <- NULL
  if (C == "pds_nP"){
    events=data.frame(var = "Fad", time = 31, value = 0.15, method = "add")
  }
  x <- x + Xs(model0, optionsOde = optionsOde, optionsSens = optionsSens, condition = C, events = events)
}

## Data ----------------------------------------------------------------------

data <- readRDS("inst/examples/example_phytoene/data_triplicates.Rds")
if(!do.fiterrors)
  data <- readRDS("inst/examples/example_phytoene/data_sigma.Rds")


## Objective Functions -------------------------------------------------------

# Initalize parameters 
outerpars <- getSymbols(do.call(c, trafoL[conditions]))

# ratios of scaling factors were determined in a previous fit to avoid that the model mis-uses the scaling factors l1 = s3/s1; l2 = s3/s2
scalings_ls <- c(logl1_new = log(1.06), logl2_new = log(0.614),logl1 = log(1.08), logl2 = log(0.8),logl1_WT = log(1.26), logl2_WT = log(0.898), logl1_dpq = log(1),logl2_dpq = log(1))
sigma_ls <- c(0.04,0.03,0.1,0.05,0.1,0.1, 0.05,0.05)
names(sigma_ls) <- names(scalings_ls)

set.seed(2)
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
pouter <- rnorm(length(prior), prior, 1); names(pouter) <- outerpars
cOuter <- constraintL2(mu = prior, sigma = 10) # general prior for all parameters -> removed for final results
cOuterSL <- constraintL2(mu = scalings_ls[names(scalings_ls) %in% outerpars], sigma = sigma_ls[names(scalings_ls) %in% outerpars],attr.name = "dataSC")

# Objective function for trust()
obj <- normL2(data, g*x*p) + cOuterSL  #+ cOuter
if(do.fiterrors)
  obj <- normL2(data, g*x*p, err) + cOuterSL #+ cOuter
