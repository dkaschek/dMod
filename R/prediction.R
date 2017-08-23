

#' Model prediction function for ODE models. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, deriv = TRUE)} returning ODE output and sensitivities.
#' @param odemodel object of class \link{odemodel}
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' Within \code{Xs()} a \code{data.frame} of additional events is generated to 
#' reset the sensitivities appropriately, depending on the event method. 
#' ATTENTION: The addional events are not dynamically recalculated. If you call the prediction
#' function with alternative events, the prediction is fine but the sensitivities can be wrong.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @return Object of class \link{prdfn}. If the function is called with parameters that
#' result from a parameter transformation (see \link{P}), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @example inst/examples/test_blocks.R
#' @export
#' @import deSolve
Xs <- function(odemodel, forcings=NULL, events=NULL, names = NULL, condition = NULL, optionsOde=list(method = "lsoda"), optionsSens=list(method = "lsodes")) {
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  myforcings <- forcings
  myevents <- events
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  
  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar
  
  #Additional events for resetting the sensitivities when events are supplied
  myevents.addon <- NULL
  if(!is.null(myevents)) {
    myevents.addon <- lapply(1:nrow(myevents), function(i) {
      newevent <- with(myevents[i, ], {
        newvar <- sensvar[senssplit.1 == var]
        newtime <- time
        newvalue <- switch(as.character(method), replace = 0, add = 0, multiply = value)
        newmethod <- method
        if(length(newvar) > 0 && method != "add") {
          data.frame(var = newvar, time = newtime, value = newvalue, method = newmethod)
        } else {
          NULL
        }
      })
      
      return(newevent)
      
    })
    myevents.addon <- do.call(rbind, myevents.addon)
  }
  
  # Names for deriv output
  sensGrid <- expand.grid(variables, c(svariables, sparameters), stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)
  
  # Update sensNames when names are set
  select <- sensGrid[, 1] %in% names
  sensNames <- paste(sensGrid[,1][select], sensGrid[,2][select], sep = ".")  
  
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    events.addon = myevents.addon,
    optionsOde = optionsOde,
    optionsSens = optionsSens
  )
  
  P2X <- function(times, pars, deriv=TRUE){
    
    
    yini <- unclass(pars)[variables]
    mypars <- unclass(pars)[parameters]
    
    events <- controls$events
    forcings <- controls$forcings
    myevents.addon <- controls$events.addon
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    names <- controls$names
    events.extended <- rbind(events, myevents.addon)
    
    # Add event time points (required by integrator) 
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))
    
    # Sort event time points
    if (!is.null(events)) events <- events[order(events$time),]
    if (!is.null(events.extended)) events.extended <- events.extended[order(events.extended$time),]
    
    
    myderivs <- NULL
    mysensitivities <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      # loadDLL(func)
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- do.call(odeC, c(list(y = unclass(yini), times = times, func = func, parms = mypars, forcings = forc, events = list(data = events)), optionsOde))
      out <- submatrix(out, cols = c("time", names))
      #out <- cbind(out, out.inputs)
      
      
    } else {
      
      # Evaluate extended model
      # loadDLL(extended)
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      
      outSens <- do.call(odeC, c(list(y = c(unclass(yini), yiniSens), times = times, func = extended, parms = mypars, 
                                      forcings = forc, 
                                      events = list(data = events.extended)), optionsSens))
      #out <- cbind(outSens[,c("time", variables)], out.inputs)
      out <- submatrix(outSens, cols = c("time", names))
      mysensitivities <- submatrix(outSens, cols = !colnames(outSens) %in% c(variables, forcnames))
      
      
      # Apply parameter transformation to the derivatives
      variables <- intersect(variables, names)
      sensLong <- matrix(outSens[,sensNames], nrow = dim(outSens)[1]*length(variables))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        sensLong <- sensLong %*% submatrix(dP, rows = c(svariables, sparameters))
        sensGrid <- expand.grid.alt(variables, colnames(dP))
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      }
      myderivs <- matrix(0, nrow = nrow(outSens), ncol = 1 + length(sensNames), dimnames = list(NULL, c("time", sensNames)))
      myderivs[, 1] <- out[, 1]
      myderivs[, -1] <- sensLong
      
    }
    
    #prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = unique(sensGrid[,2]))
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
  
  
  
  
}


#' Model prediction function for ODE models without sensitivities. 
#' @description Interface to get an ODE 
#' into a model function \code{x(times, pars, forcings, events)} returning ODE output.
#' It is a reduced version of \link{Xs}, missing the sensitivities. 
#' @param odemodel Object of class \link{odemodel}.
#' @param forcings, see \link{Xs}
#' @param events, see \link{Xs}
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to \code{f}. 
#' All quantities that are not initialised by pars 
#' in \code{x(..., forcings, events)} are initialized with 0. For more details and
#' the return value see \link{Xs}.
#' @export
Xf <- function(odemodel, forcings = NULL, events = NULL, condition = NULL, optionsOde=list(method = "lsoda")) {
  
  func <- odemodel$func
  
  myforcings <- forcings
  myevents <- events
  
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0,length(variables))
  names(yini) <- variables
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    optionsOde = optionsOde
  )
  
  P2X <- function(times, P){
    
    events <- controls$events
    forcings <- controls$forcings
    
    # Add event time points (required by integrator) 
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))
    
    
    yini[names(P[names(P) %in% variables])] <- P[names(P) %in% variables]
    pars <- P[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    # loadDLL(func)
    if(!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=pars, forcings=forc,events = list(data = events)), optionsOde))
    #out <- cbind(out, out.inputs)      
    
    prdframe(out, deriv = NULL, parameters = P)
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
  
  
  
  
  
  
  
}


#' Model prediction function from data.frame
#' 
#' @param data data.frame with columns "name", "time", and row names that 
#' are taken as parameter names. The data frame can contain a column "value"
#' to initialize the parameters.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class \link{prdfn}, i.e. 
#' a function \code{x(times pars, deriv = TRUE, conditions = NULL)}, 
#' see also \link{Xs}. Attributes are "parameters", the parameter names (row names of
#' the data frame), and possibly "pouter", a named numeric vector which is generated
#' from \code{data$value}.
#' @examples
#' # Generate a data.frame and corresponding prediction function
#' timesD <- seq(0, 2*pi, 0.5)
#' mydata <- data.frame(name = "A", time = timesD, value = sin(timesD), 
#'                      row.names = paste0("par", 1:length(timesD)))
#' x <- Xd(mydata)
#' 
#' # Evaluate the prediction function at different time points
#' times <- seq(0, 2*pi, 0.01)
#' pouter <- structure(mydata$value, names = rownames(mydata))
#' prediction <- x(times, pouter)
#' plot(prediction)
#' 
#' @export
Xd <- function(data, condition = NULL) {
  
  states <- unique(as.character(data$name))
  
  
  # List of prediction functions with sensitivities
  predL <- lapply(states, function(s) {
    subdata <- subset(data, as.character(name) == s)
    
    M <- diag(1, nrow(subdata), nrow(subdata))
    parameters.specific <- rownames(subdata)
    if(is.null(parameters.specific)) parameters.specific <- paste("par", s, 1:nrow(subdata), sep = "_")
    sensnames <- paste(s, parameters.specific, sep = ".")
    
    # return function
    out <- function(times, pars) {
      value <- approx(x = subdata$time, y = pars[parameters.specific], xout = times, rule = 2)$y
      grad <- do.call(cbind, lapply(1:nrow(subdata), function(i) {
        approx(x = subdata$time, y = M[, i], xout = times, rule = 2)$y
      }))
      colnames(grad) <- sensnames
      attr(value, "sensitivities") <- grad
      attr(value, "sensnames") <- sensnames
      return(value)
    }
    
    attr(out, "parameters") <- parameters.specific
    
    return(out)
    
  }); names(predL) <- states
  
  # Collect parameters
  parameters <- unlist(lapply(predL, function(p) attr(p, "parameters")))
  
  # Initialize parameters if available
  pouter <- NULL
  if(any(colnames(data) == "value")) 
    pouter <- structure(data$value[match(parameters, rownames(data))], names = parameters)
  
  sensGrid <- expand.grid(states, parameters, stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  
  controls <- list()  
  
  P2X <- function(times, pars, deriv=TRUE){
    
    
    predictions <- lapply(states, function(s) predL[[s]](times, pars)); names(predictions) <- states
    
    out <- cbind(times, do.call(cbind, predictions))
    colnames(out) <- c("time", states)
    
    mysensitivities <- NULL
    myderivs <- NULL
    if(deriv) {
      
      # Fill in sensitivities
      outSens <- matrix(0, nrow = length(times), ncol = length(sensNames), dimnames = list(NULL, c(sensNames)))
      for(s in states) {
        mysens <- attr(predictions[[s]], "sensitivities")
        mynames <- attr(predictions[[s]], "sensnames")
        outSens[, mynames] <- mysens
      }
      
      mysensitivities <- cbind(time = times, outSens)
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens, nrow = nrow(outSens)*length(states))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        sensLong <- sensLong %*% submatrix(dP, rows = parameters)
        sensGrid <- expand.grid.alt(states, colnames(dP))
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      }
      outSens <- cbind(times, matrix(sensLong, nrow = dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      myderivs <- outSens
      #attr(out, "deriv") <- outSens
    }
    
    #attr(out, "parameters") <- unique(sensGrid[,2])
    
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions. 
#' 
#' @description Creates an object of type \link{obsfn} that evaluates an observation function
#' and its derivatives based on the output of a model prediction function, see \link{prdfn}, 
#' as e.g. produced by \link{Xs}.
#' @param g Named character vector or equation vector defining the observation function
#' @param f Named character of equations or object that can be converted to eqnvec or object of class fn.
#' If f is provided, states and parameters are guessed from f.
#' @param states character vector, alternative definition of "states", usually the names of \code{f}. If both,
#' f and states are provided, the states argument overwrites the states derived from f.
#' @param parameters character vector, alternative definition of the "parameters",
#' usually the symbols contained in "g" and "f" except for \code{states} and the code word \code{time}. If both,
#' f and parameters are provided, the parameters argument overwrites the parameters derived from f and g.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param attach.input logical, indiating whether the original input should be
#' returned with the output.
#' @param compile Logical, compile the function (see \link{funC0})
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to R command line.
#' @return Object of class \link{obsfn}, i.e.
#' a function \code{y(..., deriv = TRUE, conditions = NULL)} representing the evaluation of the 
#' observation function. Arguments \code{out} (model prediction) and \code{pars} (parameter values)
#' shoudl be passed by the \code{...} argument.
#' If \code{out} has the attribute  "sensitivities", the result of
#' \code{y(out, pars)}, will have an attributed "deriv" which reflecs the sensitivities of 
#' the observation with respect to the parameters.
#' If \code{pars} is the result of a parameter transformation \code{p(pars)} (see \link{P}), 
#' the Jacobian 
#' of the parameter transformation and the sensitivities of the observation function
#' are multiplied according to the chain rule for differentiation.
#' @example inst/examples/prediction.R
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL, condition = NULL, attach.input = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE) {
 
  
  # Idea: 
  # If replicate scaling is undispensible and different 
  # observable names for different replicates is not an option, then
  # g could be a list of observables. For this case, the observation
  # function has to return a list of observations for each condition.
  # Not yet clear how this works with the "+" operator.
  myattach.input <- attach.input
  
  warnings <- FALSE
  modelname_deriv <- NULL
  
  if (is.null(f) && is.null(states) && is.null(parameters)) 
    stop("Not all three arguments f, states and parameters can be NULL")
  
  # Modify modelname by condition
  if (!is.null(modelname) && !is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Then add suffix(es) for derivative function
  if (!is.null(modelname)) modelname_deriv <- paste(modelname, "deriv", sep = "_")
  
  # Get potential paramters from g, forcings are treated as parameters because
  # sensitivities dx/dp with respect to forcings are zero.
  if (is.null(f)) {
    states <- union(states, "time")
    parameters <- parameters
  } else if (inherits(f, "fn")) {
    states <- union(names(attr(attr(f, "mappings")[[1]], "equations")), "time")
    parameters <- setdiff(union(getParameters(f), getSymbols(unclass(g))), states)
  } else {
    f <- as.eqnvec(f)
    if (is.null(states)) states <- union(names(f), "time")
    if (is.null(parameters)) parameters <- getSymbols(c(unclass(g), unclass(f)), exclude = c(states, "time"))
  }
  variables.deriv <- c(
    states, 
    as.vector(outer(states, c(states, parameters), paste, sep = "."))
  )
  
  # Observables defined by g
  observables <- names(g)
  
  gEval <- funC0(g, variables = states, parameters = parameters, compile = compile, modelname = modelname, 
                 verbose = verbose, convenient = FALSE, warnings = FALSE)
  
  # Character matrices of derivatives
  dxdp <- dgdx <- dgdp <- NULL
  
  if (length(states) > 0 & length(parameters) > 0) {
    dxdp <- apply(expand.grid.alt(states, c(states, parameters)), 1, paste, collapse = ".")
    dxdp <- matrix(dxdp, nrow = length(states))
  }
  if (length(states) > 0)
    dgdx <- matrix(jacobianSymb(g, states), nrow = length(g))
  if (length(parameters) > 0) {
    dgdp <- cbind(
      matrix("0", nrow = length(g), ncol = length(states)), 
      matrix(jacobianSymb(g, parameters), nrow = length(g))
    )
  }
  
  # Sensitivities of the observables
  derivs <- as.vector(sumSymb(prodSymb(dgdx, dxdp), dgdp))
  if (length(derivs) == 0) stop("Neither states nor parameters involved")
  names(derivs) <- apply(expand.grid.alt(observables, c(states, parameters)), 1, paste, collapse = ".")

  derivsEval <- funC0(derivs, variables = variables.deriv, parameters = parameters, compile = compile, modelname = modelname_deriv,
                      verbose = verbose, convenient = FALSE, warnings = FALSE)
  
  # Vector with zeros for possibly missing derivatives
  # zeros <- rep(0, length(dxdp))
  # names(zeros) <- dxdp
  # Redundant -> missing values have been implemented in funC0
  
  controls <- list(attach.input = attach.input) 
  
  X2Y <- function(out, pars) {
    
    attach.input <- controls$attach.input
    
    # Prepare list for with()
    nOut <- ncol(out)
    values <- gEval(M = out, p = pars)
    
    sensitivities.export <- NULL
    dout <- attr(out, "sensitivities")
    if (!is.null(dout)) {
      dvalues <- derivsEval(M = cbind(out, dout), p = pars)
      sensitivities.export <- cbind(time = out[, 1], dvalues)
    }
    
    
    # Parameter transformation
    dP <- attr(pars, "deriv")
    if (!is.null(dP) & !is.null(dout)) {
      
      parameters.all <- c(states, parameters)
      parameters.missing <- parameters.all[!parameters.all %in% rownames(dP)]
      
      if (length(parameters.missing) > 0 & warnings)
        warning("Parameters ", paste(parameters.missing, collapse = ", ", "are missing in the Jacobian of the parameter transformation. Zeros are introduced."))
      
      dP.full <- matrix(0, nrow = length(parameters.all), ncol = ncol(dP), dimnames = list(parameters.all, colnames(dP)))
      dP.full[intersect(rownames(dP), parameters.all),] <- dP[intersect(rownames(dP), parameters.all),]

      # Multiplication with tangent map
      sensLong <- matrix(dvalues, nrow = nrow(out)*length(observables))
      sensLong <- sensLong %*% dP.full
      dvalues <- matrix(sensLong, nrow = dim(out)[1])
      
      # Naming
      sensGrid <- expand.grid.alt(observables, colnames(dP.full))
      sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      colnames(dvalues) <- sensNames
      
    }
    
    
    
    # Format output
    values <- cbind(time = out[,"time"], values)
    if (attach.input)
      values <- cbind(values, submatrix(out, cols = -1))
    
    
    myderivs <- myparameters <- NULL
    if (!is.null(dout) && !attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues)
      if (is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    if (!is.null(dout) && attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues, submatrix(attr(out, "deriv"), cols = -1))
      if (is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    
    
    # Output 
    prdframe(prediction = values, deriv = myderivs, sensitivities = sensitivities.export, parameters = pars) 
    
    
    
  }
  
  attr(X2Y, "equations") <- g
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states") <- states
  attr(X2Y, "modelname") <- modelname
  
  obsfn(X2Y, parameters, condition)
  
}

#' Generate a prediction function that returns times
#' 
#' Function to deal with non-ODE models within the framework of dMod. See example.
#' 
#' @param condition  either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class \link{prdfn}.
#' @examples 
#' x <- Xt()
#' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"))
#' 
#' times <- seq(-1, 1, by = .05)
#' pars <- c(a = .1, b = 1)
#' 
#' plot((g*x)(times, pars))
#' @export
Xt <- function(condition = NULL) {
  
  
  
  # Controls to be modified from outside
  controls <- list()
  
  P2X <- function(times, pars, deriv=TRUE){
    
    out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))
    sens <- deriv <- out
    
    prdframe(out, deriv = deriv, sensitivities = sens, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- NULL
  attr(P2X, "equations") <- NULL
  attr(P2X, "forcings") <- NULL
  attr(P2X, "events") <- NULL
  
  
  prdfn(P2X, NULL, condition) 
  
  
  
  
}


