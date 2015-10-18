## Todo:  Basic constructor for obsfn


## The obsfn class and its constructor ----------------------------------------------------------

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
  class(X2Y) <- "obsfn"
  
  return(X2Y)
  
  
}

