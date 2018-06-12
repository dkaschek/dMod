#' IQR QSP Model
#'
#' @param input character pointing to the model, see \code{\link[IQRtools]{IQRmodel}} for details.
#' @param regression character vector containing names of regression parameters.
#' @param fix character vector denoting which parameters and initial
#' values to fix during parameter estimation.
#' @param estimate character vector denoting which parameters and initial values to estimate.
#' If \code{estimate} is provided, \code{fix} is overwritten.
#' @param ndoses the maximal number of administrations to prepare for.
#' @param ... arguments going to \code{\link{odemodel}()}.
#' @return list
#' @import IQRtools
#' @export
read.IQRmodel <- function(input, regression = NULL, fix = NULL, estimate = NULL, ndoses = 1,  ...) {
  
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
  if (aux_fileparts(input)$fileext==".xml") {
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
    print(vars__)
    factors__ <- lapply(derivs__, function(x) x[x != "0"])
    print(factors__)
    
    
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
    list.files(pattern = glob2rx(paste0("par_", model__$name, "*.c*")))
  )
  output__ <- model__$name
  .so <- .Platform$dynlib.ext
  system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", paste(files__, collapse = " "), " -o ", output__, .so), intern = TRUE)
  dyn.load(paste0(output__, .so))
  modelname(p__) <- modelname(x__) <- modelname(g__) <- modelname(g_aux__) <- output__
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
  list(g = g__,  g_aux = g_aux__, x =  x__, p = p__, pars = pars__, fixed = fixed__, model = myodemodel__)
  
  
  
  
}

#' Read IQRsysData file
#' 
#' @param data character pointing to the CSV data file
#' @param regression names of the regression parameters to be imported
#' 
#' @return A datalist object with additional attributes \code{regression} and
#' \code{dosing} (dosing parameter transformations).
#' @export
read.IQRdata <- function(data, regression = NULL) {

  # Read Data
  if (is.character(data))
    mydata0 <- IQRloadCSVdata(data)
  else
    mydata0 <- data
  
  required <- c("NAME", "TIME", "DV", "CONDITION")
  conditions <- unique(mydata0[["CONDITION"]])
  mydosing <- mydata0[mydata0[["YTYPE"]] == 0, c("TIME", "CONDITION", "ID", "AMT", "ADM", "TINF")]
  mydata <- mydata0[mydata0[["YTYPE"]] != 0, ]
  
  # Generate datalist
  data <- with(mydata, {
    out1 <- data.frame(
      name = paste0("OUTPUT", YTYPE),
      time = TIME,
      value = DV,
      sigma = NA,
      condition = CONDITION
    )
    if (!is.null(regression)) {
      out1 <- cbind(out1, mydata[, regression, drop = FALSE])
    }
    out1 <- as.datalist(out1, split.by = "condition", keep.covariates = regression)
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
    dosing <- lapply(split(mydosing, mydosing[["CONDITION"]]), function(mydosing) {
      
      mydosing <- mydosing[mydosing$ID == mydosing$ID[1], c("TIME", "AMT", "ADM", "TINF")]
      if (nrow(mydosing) > 0) {
        INPUT <- mydosing[["ADM"]]
        nDOSE <- sapply(1:length(INPUT), function(i) {
          cumsum(mydosing[["ADM"]] == INPUT[i])[mydosing[["ADM"]] == INPUT[i]][i]
        })
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
  attr(output, "regression") <- regression
  return(output)
    
}

#' Export to IQRtools imports to dMod.frame format
#' 
#' @param hypothesis hypothesis name
#' @param data result from \link{read.IQRdata()}
#' @param model result from \link{read.IQRmodel()}
#' @param errmodel an error model
#' @param transform vector with transformations ("L" (log), "N" (normal) or "G" (log-it))
#' @return A \linke{dMod.frame} object
#' @export
to_dMod.frame <- function(hypothesis = date(), data, model, errmodel = NULL, transform = NULL) {
  
  conditions <- names(data)
  regression <- attr(data, "regression")
  dosing <- attr(data, "dosing")
  
  p <- Reduce("+", lapply(conditions, function(C) {
    
    trafo <- getEquations(model[["p"]])[[1]]
    # Add missing error model parameters
    errpars <- getParameters(errmodel)
    trafo[errpars[!errpars %in% names(trafo)]] <- errpars[!errpars %in% names(trafo)]
    
    covtable <- covariates(data)
    # Plug in regression parameters
    for (r in regression) {
      trafo <- replaceSymbols(r, covtable[C, r], trafo)
    }
    # Plug in dosing parameters
    mydosing <- dosing[[C]]
    trafo[names(mydosing)] <- mydosing
    
    # Replace fixed parameters in what is left
    trafo <- replaceSymbols(names(mymodel$fixed), mymodel$fixed, trafo)
    
    # Perform parameter transformations
    symbols <- getSymbols(trafo)
    tG <- names(transform)[transform == "G"]
    tN <- names(transform)[transform == "N"]
    tL <- setdiff(symbols, c(tG, tN))
    
    trafo <- repar("x ~ exp(x)", x = tL, trafo)
    trafo <- repar("x ~ exp(x)/(1+exp(x))", x = tG, trafo)
    
    P(trafo, condition = C)
    
  }))
  
  dMod.frame(hypothesis, model[["g"]], model[["x"]], p, data, errmodel)
  
  
}
