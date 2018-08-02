#' IQR QSP Model
#'
#' @param input character pointing to the model, see \code{\link[IQRtools]{IQRmodel}} for details.
#' @param errmodel eqnvec to define the error model, 
#' e.g. \code{eqnvec(OUTPUT1 = "sigma_abs_1 * OUTPUT1 + sigma_rel_1")}. If NULL, no error model is returned.
#' @param regression character vector containing names of regression parameters.
#' @param fix character vector denoting which parameters and initial
#' values to fix during parameter estimation.
#' @param estimate character vector denoting which parameters and initial values to estimate.
#' If \code{estimate} is provided, \code{fix} is overwritten.
#' @param ndoses the maximal number of administrations to prepare for.
#' @param ... arguments going to \code{\link{odemodel}()}.
#' @return list
#' @export
read_IQRmodel <- function(input, errmodel = NULL, regression = NULL, fix = NULL, estimate = NULL, ndoses = 1, ...) {
  
  
  
  # Set solver
  solver <- list(...)$solver
  if (is.null(solver)) solver <- "deSolve"
  
  # Initialize
  fix__ <- fix
  estimate__ <- estimate
  forcings__ <- regression
  
  # Auxiliary functions
  quoted__ <- function(x) paste0("(", x, ")")
  is_number <- function(x) !is.na(suppressWarnings(as.numeric(initials.values__)))
  
  # Checking of input arguments
  if (!file.exists(input)) stop("Provided input argument does not point to a file on the filesystem.")
  
  # Attempt import of model file
  if (IQRtools:::aux_fileparts(input)$fileext==".xml") {
    model__ <- IQRtools:::importSBML_IQRmodel(input)
  } else {
    model__ <- IQRtools:::import_IQRmodel(input,FLAGtextIQRmodel=FALSE)
    IQRtools:::checkNames_IQRmodel(model__)
  }
  
  
  # Get ODEs
  states__ <- names(model__[["states"]])
  states.equations__ <- sapply(states__, function(s) model__[["states"]][[s]][["ODE"]])
  # Get observables
  observables__ <- names(model__$outputs)
  observables.equations__ <- sapply(observables__, function(s) model__[["outputs"]][[s]][["formula"]])
  
  
  
  ## Replace INPUT by corresponding symbolic forcing
  # Get inputs
  symbols__ <- getSymbols(states.equations__)
  inputs__ <- symbols__[grepl("^INPUT", symbols__)]
  events__ <- NULL
  
  # Check inputs
  if (length(inputs__) > 0) {
    
    diff__ <- setdiff(inputs__, paste0("INPUT", 1:length(inputs__)))
    if (length(diff__) > 0)
      stop("Inputs must be defind in numeric order, e.g., INPUT1, INPUT2, INPUT3, ...")
    
    # Add input states if not present
    input_exists__ <- intersect(inputs__, states__)
    input_new__ <- setdiff(inputs__, states__)
    
    states__ <- union(states__, inputs__)
    if (length(input_exists__) > 0)
      states.equations__[input_exists__] <- paste(states.equations__[input_exists__], "0", sep = " + ")
    if (length(input_new__) > 0)
      states.equations__[input_new__] <- "0"
    
    # Generate events
    derivs__ <- lapply(inputs__, function(x) jacobianSymb(states.equations__, x))
    vars__ <- lapply(derivs__, function(x) states__[x != "0"])
    factors__ <- lapply(derivs__, function(x) x[x != "0"])
    
    
    events__ <- data.frame(
      var = inputs__,
      time = c(paste("ton_INPUT", 1:length(inputs__), sep = ""),
               paste("toff_INPUT", 1:length(inputs__), sep = "")),
      value = c(paste("xon_INPUT", 1:length(inputs__), sep = ""),
                rep("0", length(inputs__))),
      method = "replace",
      stringsAsFactors = FALSE
    )
    
    # Duplicate events according to number of doses  
    events__ <- do.call(rbind, lapply(seq_len(ndoses), function(n) {
      
      d <- events__
      time.is.variable <- sapply(d[["time"]], function(x) {
        x <- type.convert(as.character(x))
        !is.numeric(x)
      })
      value.is.variable <- sapply(d[["value"]], function(x) {
        x <- type.convert(as.character(x))
        !is.numeric(x)
      })
      d[["time"]][time.is.variable] <- paste(d[["time"]][time.is.variable], n, sep = "_")
      d[["value"]][value.is.variable] <- paste(d[["value"]][value.is.variable], n, sep = "_")
      return(d)
      
    }))
    
  }
  
  ## Replace variables in equations and observables
  variables__ <- names(model__$variables)
  variables.equations__ <- sapply(variables__, function(s) model__[["variables"]][[s]][["formula"]])
  variables.equations__ <- resolveRecurrence(variables.equations__)
  # Determine those replacements which depend on time
  is.timedependent__ <- sapply(variables.equations__, function(myeq) any(states__ %in% getSymbols(myeq)))
  # Replace time-dependent variables in ODEs by model variables
  states.equations__ <- replaceSymbols(variables__[is.timedependent__],
                                       quoted__(variables.equations__[is.timedependent__]),
                                       states.equations__)
  # Replace time-dependent variables in observtions
  observables.equations__ <- replaceSymbols(variables__[is.timedependent__],
                                            quoted__(variables.equations__[is.timedependent__]),
                                            observables.equations__)
  # Deal with time-independent replacements (initial values and parameter transformations)
  initials__ <- names(model__$states)
  initials.values__ <- sapply(initials__, function(s) model__[["states"]][[s]][["IC"]])
  parameters.equations__ <- c(variables.equations__[!is.timedependent__],
                              initials.values__[!is_number(initials.values__)])
  
  
  # Add parameterization of events
  
  
  ## Genererate dMod outputs
  # Set up odemodel (determine which parameters can be fixed during computation of sensitivities)
  if (!is.null(estimate__))
    fix__ <- setdiff(getSymbols(c(states.equations__, parameters.equations__), exclude = "time"), estimate__)
  myodemodel__ <- odemodel(states.equations__, estimate = estimate__, forcings = forcings__, events = events__, modelname = model__$name, compile = FALSE, ...)
  
  # Set up prediction function
  optionsOde__ <- switch(solver, "deSolve" = list(method = "lsoda"), "Sundials" = list(method = "bdf"))
  optionsSens__ <- switch(solver, "deSolve" = list(method = "lsodes"), "Sundials" = list(method = "bdf"))
  x__ <- Xs(myodemodel__, optionsOde = optionsOde__, optionsSens = optionsSens__)
  # Determine parameters in model and observation
  parameters__ <- unique(c(getParameters(myodemodel__),
                           getSymbols(observables.equations__, exclude = c(forcings__, "time")),
                           getSymbols(as.character(events__[["time"]])),
                           getSymbols(as.character(events__[["value"]]))))
  # Set up observation function
  g__ <- Y(observables.equations__,
           states = states__,
           parameters = setdiff(parameters__, c(fix__, states__)),
           attach.input = FALSE,
           compile = FALSE, modelname = paste0("obs_", model__$name))
  # Setup up output function
  g_aux__ <- Y(variables.equations__,
               states = states__,
               attach.input = FALSE, deriv = FALSE, compile = FALSE, modelname = paste0("obsaux_", model__$name))
  # Set up error model function
  e__ <- Y(errmodel, f = g__, modelname = paste0("err_", model__$name), attach.input = FALSE)
  errmodel.parameters__ <- getSymbols(errmodel, exclude = observables__)
  parameters__ <- union(parameters__, errmodel.parameters__)
  
  
  # Set up parameter transformation function
  transformation__ <- repar("x~x", x = parameters__)
  transformation__ <- repar("x~y", x = names(parameters.equations__), y = parameters.equations__, transformation__)
  
  inputparmeters__ <- NULL
  if (length(inputs__) > 0) {
    for (inp in 1:length(inputs__)) {
      for (dos in 1:ndoses) {
        transformation__ <- replaceSymbols(
          paste0("ton_INPUT", inp, "_", dos), 
          paste0("TIME", inp, "_", dos), 
          transformation__
        )
        transformation__ <- replaceSymbols(
          paste0("toff_INPUT", inp, "_", dos), 
          paste0("(TIME", inp, "_", dos, "+ TINF", inp, "_", dos, ")"), 
          transformation__
        )
        transformation__ <- replaceSymbols(
          paste0("xon_INPUT", inp, "_", dos), 
          paste0("(AMT", inp, "_", dos, "/TINF", inp, "_", dos, ")"), 
          transformation__
        )
      }
      inputparameters__ <- outer(c("TIME", "TINF", "AMT"), 1:length(inputs__), function(x, y) paste(x, y, sep = ""))
      inputparameters__ <- outer(inputparameters__, 1:ndoses, function(x, y) paste(x, y, sep = "_"))
    }
    
    # Include input parameters into the list of fixed parameters
    fix__ <- union(fix__, setdiff(inputparameters__, estimate__))
  }
  p__ <- P(transformation__, modelname = paste0("par_", model__$name), compile = FALSE)
  # Compile and load
  files__ <- c(
    list.files(pattern = glob2rx(paste0(model__$name, "*.c*"))),
    list.files(pattern = glob2rx(paste0("obs_", model__$name, "*.c*"))),
    list.files(pattern = glob2rx(paste0("obsaux_", model__$name, "*.c*"))),
    list.files(pattern = glob2rx(paste0("err_", model__$name, "*.c*"))),
    list.files(pattern = glob2rx(paste0("par_", model__$name, "*.c*")))
  )
  output__ <- model__$name
  .so <- .Platform$dynlib.ext
  system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", paste(files__, collapse = " "), " -o ", output__, .so), intern = TRUE)
  dyn.load(paste0(output__, .so))
  modelname(p__) <- modelname(x__) <- modelname(g__) <- modelname(g_aux__) <- modelname(e__)  <- output__
  # Remove source files
  unlink(paste0("*", model__$name, "*.c"))
  unlink(paste0("*", model__$name, "*.o"))
  # Create parameter vector and initalize by zero (or 1 in case of TINF)
  parameters__ <- getParameters(p__)
  parameters.values__ <- structure(rep(0, length(parameters__)), names = parameters__)
  parameters.values__[grepl("TINF", parameters__)] <- 1
  # Overwrite parameters by model definition file
  parameters__ <- intersect(names(model__$parameters), parameters__)
  parameters.values__[parameters__] <- sapply(parameters__, function(s) model__[["parameters"]][[s]][["value"]])
  parameters.values__[initials__[is_number(initials.values__)]] <- as.numeric(initials.values__[is_number(initials.values__)])
  # Split into fixed and non-fixed
  fixed__ <- parameters.values__[intersect(getParameters(p__), fix__)]
  pars__ <- parameters.values__[setdiff(names(parameters.values__), fix__)]
  
  
  # Return
  list(g = g__,  g_aux = g_aux__, e = e__, x =  x__, p = p__, pars = pars__, fixed = fixed__, model = myodemodel__)
  
  
  
  
}

#' Read IQRsysData file
#' 
#' @param data character pointing to the CSV data file
#' @param keep names of the elements to keep, e.g. regression parameters to be imported
#' @param split.by column by which to split into conditions. If other than "CONDITION" must be contained in `keep`.
#' @param loq limit of quantification, single numeric or named numerig, e.g. `c(OUTPUT1 = -3, OUTPUT2 = -5)`.
#' @return A datalist object with additional attributes \code{keep} and
#' \code{dosing} (dosing parameter transformations).
#' @export
read_IQRdata <- function(data, keep = NULL, split.by = "CONDITION", loq = -Inf) {

  # Read Data
  if (is.character(data))
    mydata0 <- IQRtools::IQRloadCSVdata(data)
  else
    mydata0 <- data
  
  required <- c("NAME", "TIME", "DV", "CONDITION")
  conditions <- unique(mydata0[["CONDITION"]])
  mydosing <- mydata0[mydata0[["YTYPE"]] == 0, c("TIME", "CONDITION", "ID", "AMT", "ADM", "TINF", keep)]
  mydata <- mydata0[mydata0[["YTYPE"]] != 0, ]
  
  # Generate datalist
  data <- with(mydata, {
    out1 <- data.frame(
      name = paste0("OUTPUT", YTYPE),
      time = TIME,
      value = DV,
      sigma = NA,
      CONDITION = CONDITION
    )
    if (!is.null(keep)) {
      out1 <- cbind(out1, mydata[, keep, drop = FALSE])
    }
    out1 <- as.datalist(out1, split.by = split.by, keep.covariates = keep)
    return(out1)
  })
  
  dosing <- NULL
  if (nrow(mydosing) > 0) {
    
    
    # Determine dosing scheme from data
    # returns a vector (length = max(ADM), value = max number of administrations)
    scheme <- Reduce(function(x, y) {
      apply(rbind(x, y), 2, max)
    }, 
    lapply(split(mydosing, mydosing[["CONDITION"]]), function(mydosing) {
      unique.dosing <- mydosing[mydosing[["ID"]] == mydosing[["ID"]][1],]
      out <- sapply(1:max(unique.dosing[["ADM"]]), function(i) {
        length(which(unique.dosing[["ADM"]] == i))
      })
      return(out)
    }))
    trafo <- as.eqnvec(do.call(c, lapply(1:length(scheme), function(ADM) {
      nDOSE <- scheme[ADM]
      names.ton <- paste0("ton_INPUT", ADM, "_", 1:nDOSE)
      ton <- paste0("TIME", ADM, "_", 1:nDOSE)
      names.toff <- paste0("toff_INPUT", ADM, "_", 1:nDOSE)
      toff <- paste0("(TIME", ADM, "_", 1:nDOSE, " + TINF", ADM, "_", 1:nDOSE, ")")
      names.xon <- paste0("xon_INPUT", ADM, "_", 1:nDOSE)
      xon <- paste0("(AMT", ADM, "_", 1:nDOSE, "/TINF", ADM, "_", 1:nDOSE, ")")
      
      setNames(c(ton, toff, xon), c(names.ton, names.toff, names.xon))
    })))
    
    # Determine dosing parameters
    dosing <- lapply(split(mydosing, mydosing[split.by]), function(mydosing) {
      
    
      mydosing <- mydosing[mydosing$ID == mydosing$ID[1], c("TIME", "AMT", "ADM", "TINF")]
      if (nrow(mydosing) > 0) {
        INPUT <- mydosing[["ADM"]]
        nDOSE <- unlist(lapply(1:length(INPUT), function(i) {
          cumsum(mydosing[["ADM"]] == INPUT[i])[mydosing[["ADM"]] == INPUT[i]]
        }))
        AMT <- mydosing[["AMT"]]
        TIME <- mydosing[["TIME"]]
        TINF <- pmax(1e-4, mydosing[["TINF"]])
        
        for (i in 1:nrow(mydosing)) {
          trafo <- replaceSymbols(paste0("TIME", INPUT[i], "_", nDOSE[i]), TIME[i], trafo)
          trafo <- replaceSymbols(paste0("TINF", INPUT[i], "_", nDOSE[i]), TINF[i], trafo)
          trafo <- replaceSymbols(paste0("AMT", INPUT[i], "_", nDOSE[i]), AMT[i], trafo)
        }
      }
      # Missing doses are placed beyond time horizon with AMT = 0
      tmax <- max(as.data.frame(data)[["time"]])*1.5
      symbols <- getSymbols(trafo)
      timepars <- symbols[grepl("TIME", symbols)]
      amtpars <- symbols[grepl("AMT", symbols)]
      tinfpars <- symbols[grepl("TINF", symbols)]
      trafo <- replaceSymbols(timepars, tmax, trafo)
      trafo <- replaceSymbols(amtpars, 0, trafo)
      trafo <- replaceSymbols(tinfpars, 1, trafo)
      
      
      
      return(trafo)
      
      
    })
    
  }
  
  output <- data
  attr(output, "dosing") <- dosing
  attr(output, "keep") <- keep
  attr(output, "loq") <- loq
  return(output)
    
}



#' Define parameterization of IQR model
#' 
#' @param guess initial guesses for parameters, named numeric.
#' @param estimate estimate paramter (1, default) or not (0), named numeric/logical.
#' @param transform parameter transformation, no ("N"), log ("L", default), or logit ("G").
#' @param iiv parameters to be estimated per individuum, subset of the parameter names. By default none.
#' @return List with guess, estimate, transform and iiv, filled up to the correct lengths. 
#' @export
define_IQRpars <- function(guess, estimate = NULL, transform = NULL, iiv_guess = NULL, iiv_estimate = NULL) {
  
  myguess <- guess
  n <- names(myguess)
  
  myestimate <- setNames(rep(1, length(n)), n)
  mytransform <- setNames(rep("L", length(n)), n)
  myiiv_guess <- setNames(rep(1, length(n)), n)
  myiiv_estimate <- setNames(rep(0, length(n)), n)
  
  if (!is.null(estimate)) {
    myestimate[intersect(names(estimate), n)] <- estimate[intersect(names(estimate), n)]
  }
  
  if (!is.null(transform)) {
    mytransform[intersect(names(transform), n)] <- transform[intersect(names(transform), n)]
  }
  
  if (!is.null(iiv_guess)) {
    myiiv_guess[intersect(names(iiv_guess), n)] <- iiv_guess[intersect(names(iiv_guess), n)]
  }
  
  if (!is.null(iiv_estimate)) {
    myiiv_estimate[intersect(names(iiv_estimate), n)] <- iiv_estimate[intersect(names(iiv_estimate), n)]
  }
  
  
  list(guess = myguess, estimate = myestimate, 
       transform = mytransform, 
       iiv_guess = myiiv_guess, iiv_estimate = myiiv_estimate)
  
}


#' Create a dMod.frame from IQRtools imports
#' 
#' @param model result from \link{read.IQRmodel()}
#' @param data result from \link{read.IQRdata()}
#' @param parameters result from \link{define.IQRpars()}
#' @param projectPath the path where to generate the project and store the results
#' @return A \link{dMod.frame} object
#' @export
IQRsysProject <- function(model, data, parameters, projectPath = NULL) {
  
  conditions <- names(data)
  keep <- attr(data, "keep")
  dosing <- attr(data, "dosing")
  loq <- attr(data, "loq")
  
  covtable <- covariates(data)
  
  fixed <- parameters$guess[parameters$estimate == 0]
  p <- Reduce("+", lapply(conditions, function(C) {
    
    trafo <- getEquations(model[["p"]])[[1]]
    
    # Plug in keep parameters
    for (r in keep) {
      trafo <- replaceSymbols(r, covtable[C, r], trafo)
    }
    # Plug in dosing parameters
    mydosing <- dosing[[C]]
    trafo[names(mydosing)] <- mydosing
    
    # Replace fixed parameters from model in what is left
    trafo <- replaceSymbols(setdiff(names(model$fixed), names(fixed)), 
                            model$fixed[setdiff(names(model$fixed), names(fixed))], 
                            trafo)
    
    # Perform parameter transformations
    symbols <- getSymbols(trafo)
    tG <- names(parameters$transform)[parameters$transform == "G"]
    tN <- names(parameters$transform)[parameters$transform == "N"]
    tL <- setdiff(symbols, c(tG, tN))
    
    trafo <- repar("x ~ exp(x)", x = tL, trafo)
    trafo <- repar("x ~ exp(x)/(1+exp(x))", x = tG, trafo)
    
    # Insert individualized parameters
    if (any(parameters$iiv_estimate > 0)) {
      trafo <- repar("x ~ (x + eta_x_condition)", x = names(parameters$iiv_guess)[parameters$iiv_estimate > 0], condition = C, trafo)
    }
    
    P(trafo, condition = C)
    
  }))
  
  
  outerpars <- getParameters(p)
  
  etapars.est <- with(parameters, {
    x <- lapply(names(iiv_guess)[iiv_estimate == 1], function(myiiv) paste("eta", myiiv, conditions, sep = "_"))
    names(x) <- names(iiv_guess)[iiv_estimate == 1]
    return(x)
  })
  
  etapars.fix <- with(parameters, {
    x <- lapply(names(iiv_guess)[iiv_estimate == 2], function(myiiv) paste("eta", myiiv, conditions, sep = "_"))
    names(x) <- names(iiv_guess)[iiv_estimate == 2]
    return(x)
  })
  
  omegapars.est <- with(parameters, {
    lapply(names(iiv_guess)[iiv_estimate == 1], function(myiiv) {
      setNames(rep(paste0("omega_", myiiv), length(etapars.est[[myiiv]])), etapars.est[[myiiv]])
    })
  })
  
  omegapars.fix <- with(parameters, {
    lapply(names(iiv_guess)[iiv_estimate == 2], function(myiiv) {
      setNames(rep(paste0("omega_", myiiv), length(etapars.fix[[myiiv]])), etapars.fix[[myiiv]])
    })
  })
  
 
  
  etanames.est <- unlist(etapars.est, use.names = FALSE)
  etanames.fix <- unlist(etapars.fix, use.names = FALSE)
  omeganames.est <- unique(unlist(omegapars.est, use.names = FALSE))
  omeganames.fix <- unique(unlist(omegapars.fix, use.names = FALSE))
  
  
  eta.est <- setNames(rep(0, length(etanames.est)), etanames.est)
  eta.fix <- setNames(rep(0, length(etanames.fix)), etanames.fix)
  omega.est <- with(parameters, setNames(log(iiv_guess[iiv_estimate == 1]), omeganames.est))
  omega.fix <- with(parameters, setNames(iiv_guess[iiv_estimate == 2], omeganames.fix))
  
  # Other paramters
  pars <- parameters$guess
  pars[parameters$transform == "L"] <- log(pars[parameters$transform == "L"])
  pars[parameters$transform == "G"] <- IQRtools::logit(pars[parameters$transform == "G"])
  fixed <- pars[parameters$estimate == 0]
  pars <- pars[parameters$estimate == 1]
  
  # Objective functions
  prd <- 
    model$g*model$x*p
  
  nlme <- NULL
  if (length(omegapars.est) > 0) {
    nlme.est <- constraintL2(mu = eta.est, sigma = unlist(omegapars.est, use.names = FALSE), attr.name = "data")
    nlme <- nlme.est
  }
  
  if (length(omegapars.fix) > 0) {
    nlme.fix <- constraintL2(mu = eta.fix, sigma = as.numeric(omega.fix[unlist(omegapars.fix, use.names = FALSE)]), attr.name = "data")
    nlme <- nlme.fix
  }
  
  if (length(omega.est) > 0 & length(omega.fix) > 0) {
    nlme <- nlme.est + nlme.fix
  }

  obj <- obj_data <- 
    normL2(data, prd, model$e, loq = loq)
  
  if (!is.null(nlme)) {
    obj <- obj_data + nlme
  }
  
  # Get data times
  timerange <- range(as.data.frame(data)$time)
  times <- seq(min(0, timerange[1]), timerange[2], length.out = 200)
  
  # Output dmod frame
  myframe <- dMod.frame(hypothesis = date(), 
                        g = model[["g"]], 
                        x = model[["x"]], 
                        p = p, 
                        data = data, 
                        e = model[["e"]])
  
  # Append objective function
  myframe <- appendObj(myframe, 
                       prd = list(prd),
                       obj_data = list(obj_data),
                       obj = list(obj),
                       parameters = list(parameters),
                       pars = list(c(omega.est, c(pars, eta.est, eta.fix))),
                       fixed = list(fixed),
                       times = list(times),
                       eta = list(c(etapars.est, etapars.fix)),
                       omega = list(omegapars.est),
                       projectPath = list(projectPath))
  
  # Copy everything to project folder
  if (!is.null(projectPath)) {
    
    if (dir.exists(projectPath)) unlink(projectPath, recursive = TRUE)
    dir.create(projectPath, recursive = TRUE)
    saveRDS(myframe, file = file.path(projectPath, "project.rds"))
    dll <- paste0(modelname(model$x), c(".so", ".dll"))
    dll <- dll[file.exists(dll)]
    
    file.copy(dll, projectPath, overwrite = TRUE)
    file.remove(dll)
    
  }
  
  return(myframe)
  
  
}


#' Runs an IQRsysProject
#' 
#' @param proj an IQRsys Project
#' @param ncores number of cores to use (uses mclapply -> only on Unix/Mac)
#' @param opt.nfits Number of fits
#' @param opt.sd Standard deviation of inital guesses for multi-start fitting
#' @param opt.iterlim Number of iterations allows for the optimizer
#' @param opt.prior_sigma Width of the L2 prior. No prior if NULL (default)
#' @param FLAGprofileLL Computes profile likelihood if TRUE (defaults to FALSE)
#' 
#' @return dMod.frame including the original project and all the results. In addition, results are written to disc.
#' 
#' @export
run_IQRsysProject <- function(proj, ncores = 1, opt.nfits = 10, opt.sd = 1, opt.iterlim = 100, opt.prior_sigma = NULL, FLAGprofileLL = FALSE) {
  
  if (.Platform$OS.type=="windows") {
    ncores <- 1
    message("Only one core supported on Windows")
  }
    
  mywd <- getwd()
  setwd(proj[["projectPath"]][[1]])
  
  # Load project from disk
  myframe <- readRDS("project.rds")
  
  # Load all available dll's
  dlls <- c(
    list.files(pattern = glob2rx("*.so")),
    list.files(pattern = glob2rx("*.dll"))
  )
  for (d in dlls) dyn.load(d)
  
  
  # Prior
  if (!is.null(opt.prior_sigma)) {
    mu <- myframe[["parameters"]][[1]][["guess"]]
    myframe <- mutate(myframe, obj = list(obj + constraintL2(mu, sigma = opt.prior_sigma)))
  }
  
  # Run fits
  cat("Running fits ... ")
  myframe <- mutate(myframe, fits = list(mstrust(obj, pars, rinit = .1, rmax = 10, sd = opt.sd, fits = opt.nfits, iterlim = opt.iterlim, cores = ncores, fixed = fixed)))
  cat("done.\n")
  
  # Augmenting by additional derived quantities
  cat("Augmenting result by derived quantities ... ")
  myframe <- appendParframes(myframe)
  myframe <- mutate(myframe, parframes_full = list(parframes))
  myframe <- mutate(myframe, parframes = list(parframes[ , !grepl("^eta_", names(parframes))]))
  myframe <- mutate(myframe, bestfit = list(as.parvec(parframes_full)))
  myframe <- mutate(myframe, vcov = list(structure(MASS::ginv(obj(bestfit, fixed = fixed)[["hessian"]]), dimnames = list(names(bestfit), names(bestfit)))))
  cat("done.\n")
  
  
  # Run profile likelihood
  if (FLAGprofileLL) {
    cat("Running profile likelihood ... ")
    myframe <- mutate(
      myframe, 
      profiles = list(
        profile(obj, bestfit, 
                whichPar = names(parameters[["guess"]])[parameters[["estimate"]] == 1 & 
                                                          !grepl("^sigma_", names(parameters[["guess"]])) &
                                                          !grepl("^omega_", names(parameters[["guess"]]))],
                cores = ncores, verbose = FALSE, fixed = fixed
        )
      )
    )
    cat("done.\n")
    
  }
  
  # Generate table of fitted parameters with uncertainty
  bestfit <- myframe[["bestfit"]][[1]]
  sigma <- sqrt(diag(myframe[["vcov"]][[1]]))
  
  mypartable <- data.frame(
    parameter = names(bestfit),
    value = as.numeric(bestfit),
    lower.68 = as.numeric(bestfit) - sigma,
    upper.68 = as.numeric(bestfit) + sigma,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # If profiles have been computed, replace sigma by that of profile
  if (FLAGprofileLL) {
    CI <- confint(myframe[["profiles"]][[1]], level = 0.68, val.column = "data")
    select <- match(CI[["name"]], mypartable[["parameter"]])
    mypartable[["lower.68"]][select] <- CI[["lower"]]
    mypartable[["upper.68"]][select] <- CI[["upper"]]
  }
  
  # Convert parameters back to original scale
  transform <- myframe[["parameters"]][[1]][["transform"]]
  estimate <- myframe[["parameters"]][[1]][["estimate"]]
  transform <- transform[estimate == 1]
  for (i in 1:length(transform)) {
    n <- names(transform)[i]
    select <- which(mypartable[["parameter"]] == n)
    v <- unlist(mypartable[select, -1])
    replacement <- switch(transform[i], N = v, L = exp(v), G = IQRtools::inv_logit(v))
    mypartable[select, -1] <- replacement
  }
  
  # Add fixed parameters
  hypothesis <- lapply(myframe[1,], function(x) x[[1]])
  mypartable[["estimate"]] <- TRUE
  fixedtable <- with(hypothesis, {
    isFixed <- parameters$estimate == 0
    data.frame(
      parameter = names(parameters$guess)[isFixed],
      value = (parameters$guess)[isFixed],
      lower.68 = NaN,
      upper.68 = NaN,
      estimate = FALSE,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  
  mypartable <- rbind(mypartable, fixedtable)
  
  # Sort table alphabetically
  mypartable <- mypartable[order(mypartable[["parameter"]]),]
  
  
  myframe <- mutate(myframe, partable = list(mypartable))
  myframe <- mutate(myframe, partable_pop = list(with(mypartable, mypartable[!grepl("^omega_", parameter) & !grepl("^eta_", parameter) & !grepl("^sigma_", parameter),])))
  myframe <- mutate(myframe, partable_err = list(with(mypartable, mypartable[grepl("^sigma_", parameter),])))
  myframe <- mutate(myframe, partable_eta = list(with(mypartable, mypartable[grepl("^eta_", parameter),])))
  myframe <- mutate(myframe, partable_omega = list(with(mypartable, mypartable[grepl("^omega_", parameter),])))
  
  # Produce tables
  cat("Producing tables ... ")
  mytables <- c("partable_pop", "partable_eta", "partable_omega", "partable_err")
  for (tabname in mytables) {
    tab <- myframe[[tabname]][[1]]
    IQRtools::IQRexportCSVdata(tab, filename = file.path("RESULTS", paste0(tabname, ".csv")))
    IQRtools::IQRoutputTable(tab, filename = file.path("RESULTS", paste0(tabname, ".txt")), report = FALSE)
  }
  cat("done.\n")
  
  # Produce plots
  cat("Producing plots ... ")
  
  # Waterfall plot
  p <- plotValues(myframe[["parframes"]][[1]]) +
    xlab("index of fit sorted by objective value") +
    ylab("objective value") +
    ggtitle("Objective values after parameter estimation for multiple initial guesses")
  IQRtools::IQRoutputPDF(p, file.path("RESULTS", "plotValues.pdf"))
  
  # Parameter values
  p <- plotPars(myframe[["parframes"]][[1]]) + xlab(NULL) + ylab("value") + coord_flip() + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ggtitle("Values of the estimated parameters for multiple initial guesses")
  IQRtools::IQRoutputPDF(p, file.path("RESULTS", "plotPars.pdf"))
  
  # Model prediction vs data
  hypothesis <- lapply(myframe[1,], function(x) x[[1]])
  p <- suppressMessages(with(hypothesis, {
    prediction <- as.data.frame(prd(times, bestfit, fixed = fixed), data = data, errfn = e) 
    
    IDs <- unique(prediction[["ID"]])
    IDs <- split(IDs, ceiling(seq_along(IDs)/12))
    
    lapply(IDs, function(myID) {
      myprediction <- prediction[prediction[["ID"]] %in% myID, ]
      mydata <- as.data.frame(data)
      mydata <- mydata[mydata[["ID"]] %in% myID, ]

      ggplot(myprediction, aes(x = time, y = value)) + 
        facet_wrap(~ID, scales = "free", nrow = 3, ncol = 4) +
        geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = .3, lty = 0) +
        geom_line() +
        geom_point(data = mydata) +
        theme_dMod()
      
    })
  }))
  IQRtools::IQRoutputPDFstart(file.path("RESULTS", "plotPredVsData.pdf"))
  for (i in 1:length(p)) print(p[[i]])
  IQRtools::IQRoutputPDFend(file.path("RESULTS", "plotPredVsData.pdf"))
  
  # Profile likelihood
  if (FLAGprofileLL) {
    p <- plotProfile(myframe[["profiles"]][[1]], mode == "data") +
      theme(legend.position = "none") +
      ggtitle("Profile likelihood")
    IQRtools::IQRoutputPDF(p, file.path("RESULTS", "plotProfile.pdf"))
  }
  
  cat("done.\n")
  
  # Save results
  saveRDS(myframe, file = "project.rds")
  
  # Return to original working directory
  setwd(mywd)
  
  # Return the results as dMod.frame
  return(myframe)
  
  
}
