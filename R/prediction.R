
#' Model evaluation. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, forcings, events, deriv = TRUE)} returning ODE output and sensitivities.
#' @param func return value from \code{funC(f)} where \code{f} defines the ODE.
#' @param extended return value from \code{funC(c(f, sensitivitiesSymb(f)))}.
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' Within \code{Xs()} a \code{data.frame} of additional events is generated to 
#' reset the sensitivities appropriately, depending on the event method. 
#' ATTENTION: The addional events are not dynamically recalculated. If you call the prediction
#' function with alternative events, the prediction is fine but the sensitivities can be wrong.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @return A model prediction function \code{x(times, pars, forcings, events, deriv = TRUE)} representing 
#' the model evaluation. The result of
#' \code{x(times, pars, forcings, events, deriv = TRUE)} contains
#' attributes "sensitivities" and "deriv" with the sensitivities if \code{deriv=TRUE}. 
#' If \code{deriv=FALSE}, sensitivities are not computed (saving time).
#' If \code{pars} is
#' the result of \code{p(pouter)} (see \link{P}), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @export
Xs <- function(func, extended, forcings=NULL, events=NULL, optionsOde=list(method="lsoda"), optionsSens=list(method="lsodes")) {
  
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
  
  P2X <- function(times, pars, forcings = NULL, events = NULL, deriv=TRUE){
    
    yini <- pars[variables]
    mypars <- pars[parameters]
    
    if(is.null(forcings)) forcings <- myforcings
    if(is.null(events)) events <- myevents
    
    myderivs <- NULL
    mysensitivities <- NULL
    if(!deriv) {
      
      # Evaluate model without sensitivities
      loadDLL(func)
      if(!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc, events = list(data = events)), optionsOde))
      #out <- cbind(out, out.inputs)
      
      
    } else {
      
      # Evaluate extended model
      loadDLL(extended)
      if(!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      outSens <- do.call(odeC, c(list(y=c(yini, yiniSens), times=times, func=extended, parms=mypars, 
                                      forcings=forc, 
                                      events = list(data = rbind(events, myevents.addon))), optionsSens))
      #out <- cbind(outSens[,c("time", variables)], out.inputs)
      out <- outSens[,c("time", c(variables, forcnames))]
      mysensitivities <- submatrix(outSens, cols = !colnames(outSens)%in%c(variables, forcnames))
      
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens[,sensNames], nrow=dim(outSens)[1]*length(variables))
      dP <- attr(pars, "deriv")
      if(!is.null(dP)) {
        sensLong <- sensLong%*%submatrix(dP, rows = c(svariables, sparameters))
        sensGrid <- expand.grid(variables, colnames(dP), stringsAsFactors=FALSE)
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      }
      outSens <- cbind(outSens[,1], matrix(sensLong, nrow=dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      myderivs <- outSens
      
    }
   
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = unique(sensGrid[,2]))
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  
  return(P2X)
  
  
}


#' Model evaluation without sensitivities. 
#' @description Interface to get an ODE 
#' into a model function \code{x(times, pars, forcings, events)} returning ODE output.
#' It is a reduced version of \link{Xs}, missing the sensitivities. 
#' @param func return value from \code{funC(f)} where \code{f} defines the ODE. 
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to \code{f}. All quantities that are not initialised by pars 
#' in \code{x(times, pars, forcings, events)} are initialized at 0.
Xf <- function(func, forcings=NULL, events=NULL, optionsOde=list(method="lsoda")) {
  
  myforcings <- forcings
  myevents <- events
  
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0,length(variables))
  names(yini) <- variables
  
  P2X <- function(times, P, forcings = myforcings, events = myevents){
    
    yini[names(P[names(P) %in% variables])] <- P[names(P) %in% variables]
    pars <- P[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    loadDLL(func)
    if(!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=pars, forcings=forc,events = list(data = myevents)), optionsOde))
    #out <- cbind(out, out.inputs)      
      
    prdframe(out, deriv = NULL, parameters = names(P))
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  return(P2X)
  
}


#' Model prediction function from data.frame
#' 
#' @param data data.frame with columns "name", "times", and row names that 
#' are taken as parameter names. The data frame can contain a column "value"
#' to initialize the parameters.
#' @return Prediction function, a function \code{x(times pars, deriv = TRUE)}, 
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
#' pouter <- attr(x, "pouter")
#' prediction <- x(times, pouter)
#' plotPrediction(list(prediction = prediction))
#' 
#' # Evaluate the sensitivities at these time points
#' sensitivities <- attr(prediction, "deriv")
#' plotPrediction(list(sens = sensitivities))
#' 
#' @export
#' 
Xd <- function(data) {
  
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
      if(!is.null(dP)) {
        sensLong <- sensLong%*%submatrix(dP, rows = parameters)
        sensGrid <- expand.grid(states, colnames(dP), stringsAsFactors = FALSE)
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      }
      outSens <- cbind(times, matrix(sensLong, nrow=dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      myderivs <- outSens
      #attr(out, "deriv") <- outSens
    }
   
    #attr(out, "parameters") <- unique(sensGrid[,2])
    
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = unique(sensGrid[,2]))
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  return(P2X)
  
}


#' Observation functions. 
#' @description Creates a function \code{y(out, pars)} that evaluates an observation function
#' and its derivatives based on the output of a model function \code{x(times, pars)}, see \link{Xf} and \link{Xs}.
#' @param g Named character vector defining the observation function
#' @param f Named character, the underlying ODE
#' @param compile Logical, compile the function (see \link{funC0})
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @return a function \code{y(out, pars, attach.input = FALSE)} representing the evaluation of the 
#' observation function. 
#' If \code{out} has the attribute  "sensitivities", the result of
#' \code{y(out, pars)}, will have an attributed "deriv" which reflec the sensitivities of 
#' the observation with respect to the parameters.
#' If \code{pars} is the result of a parameter transformation \code{p(pars)} (see \link{P}), 
#' the Jacobian 
#' of the parameter transformation and the sensitivities of the observation function
#' are multiplied according to the chain rule for differentiation.
#' If \code{attach.input = TRUE}, the original argument \code{out} will be attached to the evaluated observations.
#' @export
Y <- function(g, f, states = NULL, parameters = NULL, compile = FALSE, modelname = NULL) {
  
  warnings <- FALSE
  modelname_deriv <- NULL
  if(!is.null(modelname)) modelname_deriv <- paste(modelname, "deriv", sep = "_")
  
  # Get potential paramters from g, forcings are treated as parameters because
  # sensitivities dx/dp with respect to forcings are zero.
  if(is.null(f)) {
    states <- states
    parameters <- parameters
  } else {
    states <- names(f)
    parameters <- getSymbols(c(g, f), exclude = c(states, "time"))
  }
  
  # Observables defined by g
  observables <- names(g)
  
  gEval <- funC0(g, compile = compile, modelname = modelname)
  
  # Character matrices of derivatives
  dxdp <- dgdx <- dgdp <- NULL
  
  if(length(states) > 0 & length(parameters) > 0) {
    dxdp <- apply(expand.grid.alt(states, c(states, parameters)), 1, paste, collapse = ".")
    dxdp <- matrix(dxdp, nrow = length(states))
  }
  if(length(states) > 0)
    dgdx <- matrix(jacobianSymb(g, states), nrow=length(g))
  if(length(parameters) > 0) {
    dgdp <- cbind(
      matrix("0", nrow=length(g), ncol=length(states)), 
      matrix(jacobianSymb(g, parameters), nrow=length(g))
    )
  }
  
  # Sensitivities of the observables
  derivs <- as.vector(sumSymb(prodSymb(dgdx, dxdp), dgdp))
  if(length(derivs) == 0) stop("Nor states or parameters involved")
  names(derivs) <- apply(expand.grid.alt(observables, c(states, parameters)), 1, paste, collapse = ".")
  
  
  derivsEval <- funC0(derivs, compile = compile, modelname = modelname_deriv)
  
  # Vector with zeros for possibly missing derivatives
  zeros <- rep(0, length(dxdp))
  names(zeros) <- dxdp
  
  
  X2Y <- function(out, pars, attach.input = FALSE) {
    
    
    
    # Prepare list for with()
    nOut <- dim(out)[2]
    outlist <- lapply(1:nOut, function(i) out[,i]); names(outlist) <- colnames(out)
    
    dout <- attr(out, "sensitivities")
    if(!is.null(dout)) {
      nDeriv <- dim(dout)[2]
      derivlist <- lapply(1:nDeriv, function(i) dout[,i]); names(derivlist) <- colnames(dout)  
    } else {
      derivlist <- NULL
    }
    
    
    x <- c(outlist, derivlist, as.list(pars), as.list(zeros))
    
    values <- gEval(x)
    if(!is.null(dout)) dvalues <- derivsEval(x)
    
    # Parameter transformation
    dP <- attr(pars, "deriv")
    if(!is.null(dP) & !is.null(dout)) {
      
      parameters.all <- c(states, parameters)
      parameters.missing <- parameters.all[!parameters.all%in%rownames(dP)]
      
      if(length(parameters.missing) > 0 & warnings)
        warning("Parameters ", paste(parameters.missing, collapse = ", ", "are missing in the Jacobian of the parameter transformation. Zeros are introduced."))
      
      dP.missing <- matrix(0, nrow = length(parameters.missing), ncol=dim(dP)[2], 
                           dimnames=list(parameters.missing, colnames(dP)))
      dP <- rbind(dP, dP.missing)
      
      # Multiplication with tangent map
      sensLong <- matrix(dvalues, nrow=dim(out)[1]*length(observables))
      sensLong <- sensLong%*%submatrix(dP, rows = parameters.all)
      dvalues <- matrix(sensLong, nrow=dim(out)[1])
      
      # Naming
      sensGrid <- expand.grid.alt(observables, colnames(dP))
      sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      colnames(dvalues) <- sensNames
      
    }
    
    
    
    # Format output
    values <- cbind(time = out[,"time"], values)
    if(attach.input)
      values <- cbind(values, submatrix(out, cols = -1))
    
    
    myderivs <- myparameters <- NULL
    if(!is.null(dout) && !attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues)
      if(is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    if(!is.null(dout) && attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues, submatrix(attr(out, "deriv"), cols = -1))
      if(is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    
    
    # Output 
    prdframe(prediction = values, deriv = myderivs, parameters = myparameters) 
    
    
    
  }
  
  attr(X2Y, "equations") <- g
  
  return(X2Y)
  
  
}



#' Model evaluation. 
#' @description Interface to combine an ODE, its adjoint sensitivity equations as well
#' as the chi² and the gradient of chi²
#' into one model function \code{x(times, pars, ...)}.
#' @param func_fa Return value of \code{funC(fa)} where \code{fa} defines the ODE and the 
#' adjoint sensitivities. (See also \link{generateModelIE}).
#' @param func_l Return value of \code{funC(l))} where \code{l} defines the ODE for the
#' time-continuous chi² function and its gradient. (See also \link{generateModelIE}).
#' @param forcings A data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events (Not yet functional in bvptwp).
#' @param data A data.frame with columns "name", "time", "value" and "sigma".
#' @param optionsBvp list with arguments to be passed to \link{bvptwpC} for solving the BVP problem func_fa.
#' @param optionsOde list with arguments to be passed to \link{odeC} for the ODE integration of func_l.
#' @param optionsData list with arguments to be passed to \link{data2forc} for generation of
#' the data representation functions.
#' @return a function \code{x(times, pars, forcings, events, eps, atol, nmax, guess, deriv)} representing the model evaluation. 
#' The default values are \code{eps = 1}, \code{atol = 1e-4}, \code{nmax = 50*length(times)}, \code{guess = NULL}
#' and \code{deriv = TRUE}.
#' The result of
#' \code{x(times, pars, ...)} has attributes "value" (the chi² value) 
#' and "grad" (gradient of the chi² value) when \code{deriv=TRUE}. 
#' If \code{deriv=FALSE}, the time-continuous chi² is not evaluated.
#' If \code{pars} is 
#' the result of \code{p(P)} (see \link{P}), the Jacobian of the parameter transformation
#' and the gradient of the chi² are multiplied according to the chain rule for
#' differentiation.
Xv <- function(func_fa, func_l, forcings=NULL, events=NULL, data, optionsBvp=NULL, optionsOde = list(method="lsode"), optionsData = list()) {
  
  myforcings <- forcings
  myevents <- events
  
  variables <- attr(func_fa, "variables")
  isKnown <- which(variables%in%data$name)
  
  variablesD <- paste0(variables, "D")
  variablesD <- variablesD[isKnown]
  parameters <- attr(func_fa, "parameters")
  boundary <- attr(func_fa, "boundary")
  u <- attr(func_fa, "inputs")
  
  variables_l <- attr(func_l, "variables")
  
  dataforcings <- do.call(data2forc, c(list(data = data), optionsData))
  
  
  P2X <- function(times, pars, forcings = myforcings, events = myevents, eps = 1, atol=1e-4, nmax = 50*length(times), guess=NULL, deriv = TRUE){
    
    # Initials and parameters
    yini <- yend <- mypars <- NULL
    
    if(length(variables%in%names(pars))>0) {
      yiniorend <- pars[variables[variables%in%names(pars)]]
      yini <- yiniorend[!is.na(boundary$yini[match(names(yiniorend), boundary$name)])]
      #yend <- yiniorend[!is.na(boundary$yend[match(names(yiniorend), boundary$name)])]
    }
    if(length(pars[parameters]>0)) mypars <- pars[parameters]; mypars <- mypars[!is.na(mypars)]
    
    
    # Forcings
    forcings <- rbind(
      dataforcings, # measured time courses
      myforcings # model forcings
    )
    forc <- setForcings(func_fa, forcings)
    
    # Construct guesses 
    if(!is.null(guess)) times <- guess[,"time"]
    xguess <- times
    yguess <- matrix(0, length(variables), length(times)); rownames(yguess) <- variables
    
    if(!is.null(guess)) {
      isProvided <- variables[variables%in%colnames(guess)]
      yguess[isProvided,] <- t(guess[,isProvided])
    } else {
      guessData <- do.call(rbind, lapply(forc[variablesD], function(f) spline(f[,1], f[,2], xout=times)$y))
      yguess[isKnown,] <- guessData 
    }
    
    
    yend <- yguess[!is.na(boundary$yend[match(variables, boundary$name)]), length(xguess)]
    
    #matplot(xguess, t(yguess), type="l", lty=1)
    
    #print(dim(yguess))
    
    # Solve BVP
    
    outAll <- do.call(bvptwpC, c(list(yini = yini, 
                                      x = times, 
                                      func = func_fa, 
                                      yend = yend, 
                                      parms = c(mypars, eps = eps), 
                                      xguess = xguess, 
                                      yguess = yguess, 
                                      forcings = forc,
                                      atol = atol,
                                      allpoints = FALSE,
                                      nmax = nmax), 
                                 optionsBvp))
    colnames(outAll)[1] <- "time"
    
    # attach inputs computed from outAll
    if(!is.null(u)) {
      forcvalues <- lapply(forc, function(myforc) approxfun(myforc[,1], myforc[,2])(outAll[,1]))
      outvalues <- as.list(as.data.frame(outAll))
      u <- as.data.frame(lapply(u, function(ui) with(c(as.list(mypars), outvalues, forcvalues), eval(parse(text = ui)))))
      outAll <- cbind(outAll, u)  
    }
    
    outAll <- as.matrix(outAll)
    
    attr(outAll, "forcings") <- forc
    
    if(deriv) {
      
      ## Log-likelihood and gradient equations
      
      # Initials, parameters from above
      yini <- rep(0, length(variables_l)); names(yini) <- variables_l
      
      # Forcings
      forcings <- rbind(
        wide2long(outAll), # solution of adjoint equations
        forcings # forcings from above
      )
      forc <- setForcings(func_l, forcings)
      
      # Solve equations
      obj <- do.call(odeC, c(list(y=yini, times=outAll[,"time"], func=func_l, parms=mypars, forcings=forc), optionsOde))
      last <- dim(obj)[1]
      width <- dim(obj)[2]
      
      
      
      # Extract value and gradient
      value <- as.numeric(obj[last, 2])
      grad_par <- NULL
      grad_ini <- NULL
      
      parameters <- names(pars)[names(pars)%in%parameters]
      
      if(length(parameters)>0) {
        grad_par_names <- paste("chi", parameters, sep=".")
        grad_par <- as.numeric(obj[last, grad_par_names]); names(grad_par) <- parameters
      }
      variables <- names(pars)[names(pars)%in%variables]
      if(length(variables)>0) {
        grad_ini_names <- paste0("adj", variables)
        grad_ini <- as.numeric(2*outAll[1, grad_ini_names]); names(grad_ini) <- variables  
      }
      gradient <- c(grad_ini, grad_par)
      
      
      # Do parameter transformation
      dP <- attr(pars, "deriv")      
      
      if(!is.null(dP)) {
        gradient <- as.vector(gradient%*%(dP[names(gradient),]))
        names(gradient) <- colnames(dP)
      }
      
      attr(outAll, "value") <- value
      attr(outAll, "gradient") <- gradient
      
    }
    
    
    return(outAll)
    
  }
  
  return(P2X)
  
  
}


