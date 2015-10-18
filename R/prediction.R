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


