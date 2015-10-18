## Todo:  Decide on list of attributes that each prdfn should carry

## Class prdfn and its constructors ----------------------------------------------------


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
    myevents.addon <- lapply(1:nrow(events), function(i) {
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
  
  P2X <- function(times, pars, forcings = myforcings, events = myevents, deriv=TRUE){
    
    myforcings <- forcings
    myevents <- events
    
    yini <- pars[variables]
    mypars <- pars[parameters]
    
    myderivs <- NULL
    mysensitivities <- NULL
    if(!deriv) {
      
      # Evaluate model without sensitivities
      loadDLL(func)
      if(!is.null(myforcings)) forc <- setForcings(func, myforcings) else forc <- NULL
      out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc, events = list(data = events)), optionsOde))
      #out <- cbind(out, out.inputs)
      
      
    } else {
      
      # Evaluate extended model
      loadDLL(extended)
      if(!is.null(myforcings)) forc <- setForcings(extended, myforcings) else forc <- NULL
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
  
  class(P2X) <- "prdfn"
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
  
  P2X <- function(times, P, changedForcings = NULL, events = myevents){
    
    if(!is.null(changedForcings)) myforcings <- changedForcings
    yini[names(P[names(P) %in% variables])] <- P[names(P) %in% variables]
    pars <- P[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    loadDLL(func)
    if(!is.null(myforcings)) forc <- setForcings(func, myforcings) else forc <- NULL
    out <- do.call(odeC, c(list(y=yini, times=times, func=func, parms=pars, forcings=forc,events = list(data = events)), optionsOde))
    #out <- cbind(out, out.inputs)      
      
    prdframe(out, deriv = NULL, parameters = names(P))
    
  }
  
  class(P2X) <- "prdfn"
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
  class(P2X) <- "prdfn"
  return(P2X)
  
}

## Methods for class prdfn -------------------------------------------------



## Methods for class prdlist ------------------------------------------------

#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  prdlist(newlist)
  
}
