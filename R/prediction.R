#' Model evaluation. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, forcings, events, deriv = TRUE)} returning ODE output and sensitivities.
#' @param func return value from \code{funC(f)} where \code{f} defines the ODE.
#' @param extended return value from \code{funC(c(f, sensitivitiesSymb(f)))}.
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link{events}.
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
Xs <- function(func, extended, forcings=NULL, events=NULL, optionsOde=list(method="lsoda"), optionsSens=list(method="lsodes")) {
  
  myforcings <- forcings
  myevents <- events
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  
  
  ## Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar

  # Names for deriv output
  sensGrid <- expand.grid(variables, c(svariables, sparameters), stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  P2X <- function(times, pars, forcings = myforcings, events = myevents, deriv=TRUE){
    
    myforcings <- forcings
    myevents <- events
    
    yini <- pars[variables]
    mypars <- pars[parameters]
    
    # Interpolate forcings for output with the prediction
    out.inputs <- NULL
    if(!is.null(myforcings)) {
      inputs <- unique(myforcings$name)
      out.inputs <- unlist(lapply(inputs, function(inp) {
        t <- myforcings[myforcings$name == inp, "time"]
        y <- myforcings[myforcings$name == inp, "value"]
        approx(t, y, times)$y
      }))
      out.inputs <- matrix(out.inputs, ncol=length(inputs), dimnames = list(NULL, inputs))    
    }
    
    if(!deriv) {
    
      # Evaluate model without sensitivities
      if(!is.null(myforcings)) forc <- setForcings(func, myforcings) else forc <- NULL
      out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc, events = list(data = events)), optionsOde))
      out <- cbind(out, out.inputs)
      
      
    } else {
      
      # Evaluate extended model
      if(!is.null(myforcings)) forc <- setForcings(extended, myforcings) else forc <- NULL
      outSens <- do.call(odeC, c(list(y=c(yini, yiniSens), times=times, func=extended, parms=mypars, forcings=forc, events = list(data = events)), optionsSens))
      out <- cbind(outSens[,c("time", variables)], out.inputs)
      attr(out, "sensitivities") <- outSens[,!colnames(outSens)%in%variables]
      
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens[,sensNames], nrow=dim(outSens)[1]*length(variables))
      dP <- attr(pars, "deriv")
      if(!is.null(dP)) {
        sensLong <- sensLong%*%(dP[c(svariables, sparameters),])
        sensGrid <- expand.grid(variables, colnames(dP), stringsAsFactors=FALSE)
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      }
      outSens <- cbind(outSens[,1], matrix(sensLong, nrow=dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      
      attr(out, "deriv") <- outSens
      attr(out, "parameters") <- unique(sensGrid[,2])
    }
    
    return(out)
    
  }
  
  return(P2X)
  
  
}



#' Observation functions. 
#' @description Creates a function \code{y(out, pars)} that evaluates an observation function
#' and its derivatives based on the output of a model function \code{x(times, pars)}, see \link{X} and \link{Xs}.
#' @param g Named character vector defining the observation function
#' @param f Named character, the underlying ODE
#' @return a function \code{y(out, pars, attach=FALSE)} representing the evaluation of the observation function. 
#' If \code{out} has the attribute  "sensitivities", the result of
#' \code{y(out, pars)}, will have an attributed "deriv" which reflec the sensitivities of 
#' the observation with respect to the parameters.
#' If \code{pars} is the result of a parameter transformation \code{p(pars)} (see \link{P}), 
#' the Jacobian 
#' of the parameter transformation and the sensitivities of the observation function
#' are multiplied according to the chain rule for differentiation.
#' If \code{attach = TRUE}, the original argument \code{out} will be attached to the evaluated observations.
Y <- function(g, f) {
  
  # Get potential paramters from g, forcings are treated as parameters because
  # sensitivities dx/dp with respect to forcings are zero.
  states <- names(f)
  parameters <- getSymbols(c(g, f), exclude = c(states, "time"))
    
  # Observables defined by g
  observables <- names(g)
  gEval <- funC.algebraic(g)
  
  # Character matrices of derivatives  
  dxdp <- apply(expand.grid.alt(states, c(states, parameters)), 1, paste, collapse = ".")
  dxdp <- matrix(dxdp, nrow = length(states))
  dgdx <- matrix(jacobianSymb(g, states), nrow=length(g))
  dgdp <- cbind(
    matrix("0", nrow=length(g), ncol=length(states)), 
    matrix(jacobianSymb(g, parameters), nrow=length(g))
    )
  
  # Sensitivities of the observables
  derivs <- as.vector(sumSymb(prodSymb(dgdx, dxdp), dgdp))
  names(derivs) <- apply(expand.grid.alt(observables, c(states, parameters)), 1, paste, collapse = ".")
  derivsEval <- funC.algebraic(derivs)
    
  # Vector with zeros for possibly missing derivatives
  zeros <- rep(0, length(dxdp))
  names(zeros) <- dxdp
  
  
  X2Y <- function(out, pars, attach=FALSE) {
    
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
      
      if(length(parameters.missing) > 0)
        warning("Parameters ", paste(parameters.missing, collapse = ", ", "are missing in the Jacobian of the parameter transformation"))
      
      dP.missing <- matrix(0, nrow = length(parameters.missing), ncol=dim(dP)[2], 
                           dimnames=list(parameters.missing, colnames(dP)))
      dP <- rbind(dP, dP.missing)
      
      # Multiplication with tangent map
      sensLong <- matrix(dvalues, nrow=dim(out)[1]*length(observables))
      sensLong <- sensLong%*%(dP[parameters.all,])
      dvalues <- matrix(sensLong, nrow=dim(out)[1])
      
      # Naming
      sensGrid <- expand.grid.alt(observables, colnames(dP))
      sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")
      colnames(dvalues) <- sensNames
      
    }
  
  
  
    # Format output
    values <- cbind(time = out[,"time"], values)
    if(attach) 
      values <- cbind(values, out[,-1])
    
    
    if(!is.null(dout) && !attach) {
      attr(values, "deriv") <- cbind(time = out[,"time"], dvalues)
      if(is.null(dP)) attr(values, "parameters") <- names(pars) else attr(values, "parameters") <- colnames(dP)
    }
    if(!is.null(dout) && attach) {
      attr(values, "deriv") <- cbind(time = out[,"time"], dvalues, attr(out, "deriv")[,-1])
      if(is.null(dP)) attr(values, "parameters") <- names(pars) else attr(values, "parameters") <- colnames(dP)
    }
    return(values)        
  

  }


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
    
    myforcings <- forcings
    myevents <- events
    
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
