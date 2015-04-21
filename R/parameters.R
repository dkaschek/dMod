#' Parameter transformation
#' 
#' @param P named character vector. Names correspond to the parameter names to be fed into
#' the model. The values of P are equations that express the parameters in terms of other
#' parameters
#' @param parameters character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(P)} the identity transformation is assumed.
#' @return a function \code{p(P)} representing the parameter transformation. The result of p(P) 
#' contains an attribute "deriv" which is the jacobian of the parameter transformation
#' evaluated at \code{P}.
# P <- function(P, parameters=NULL) {
#   
#   # get outer parameters
#   fparse <- getParseData(parse(text=P))
#   symbols <- unique(fparse$text[fparse$token == "SYMBOL"])
#   if(is.null(parameters)) {
#     parameters <- symbols 
#   } else {
#     identity <- parameters[which(!parameters%in%symbols)]
#     names(identity) <- identity
#     P <- c(P, identity)
#   }
#   
#   # expresion list for parameter and jacobian evaluation
#   expressionList <- lapply(P, function(myrel) parse(text=as.character(myrel)))
#   listJac <- unlist(lapply(parameters, function(var) {
#     unlist(lapply(expressionList, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
#   }))
#   
#   jacNames <- expand.grid.alt(names(P), parameters)
#   jacNames <- paste(jacNames[,1], jacNames[,2], sep=".")
#   
#   dP <- listJac; names(dP) <- jacNames
# 
#   PEval <- funC.algebraic(P)
#   dPEval <- funC.algebraic(dP)
#   
#   # the parameter transformation function to be returned
#   p2p <- function(p, fixed=NULL, derivs = TRUE) {
#     
#     # replace fixed parameters by values given in fixed
#     if(!is.null(fixed)) {
#       is.fixed <- which(names(p)%in%names(fixed)) 
#       if(length(is.fixed)>0) p <- p[-is.fixed]
#       p <- c(p, fixed)
#     }
#     
#     # check for parameters which are not defined in parameters
#     emptypars <- names(p)[!names(p)%in%parameters]
#     
#     x <- as.list(p)
#     values <- PEval(x)[1,]
#     jacValues <- dPEval(x)
#       
#     
#     jacobian <- matrix(jacValues, ncol=length(parameters), nrow=length(P))
#     if(!is.null(fixed)) {
#       jacobian <- matrix(jacobian[,-which(parameters %in% names(fixed))], nrow=length(P))
#       colnames(jacobian) <- parameters[!parameters%in%names(fixed)]
#       rownames(jacobian) <- names(P)
#     } else {
#       colnames(jacobian) <- parameters
#       rownames(jacobian) <- names(P)
#     }
#     
#     # Append zeros for emptypars
#     if(length(emptypars) > 0) {
#       empty <- matrix(0, nrow=dim(jacobian)[1], ncol=length(emptypars))
#       colnames(empty) <- emptypars
#       jacobian <- cbind(jacobian, empty)
#     }    
#     if(derivs) attr(values, "deriv") <- jacobian
#     
#     
#     
#     return(values)
#   }
#   
#   
#   return(p2p)
#   
# }


#' Parameter transformation
#' 
#' @param trafo Named character vector. Names correspond to the parameters being fed into
#' the model (the inner parameters). The elements of tafo are equations that express 
#' the inner parameters in terms of other parameters (the outer parameters)
#' @param parameters Character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(trafo)} the identity transformation is assumed.
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @seealso \link{Pi} for implicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
P <- function(trafo, parameters=NULL) {
  
  # get outer parameters
  fparse <- getParseData(parse(text=trafo))
  symbols <- unique(fparse$text[fparse$token == "SYMBOL"])
  if(is.null(parameters)) {
    parameters <- symbols 
  } else {
    identity <- parameters[which(!parameters%in%symbols)]
    names(identity) <- identity
    trafo <- c(trafo, identity)
  }
  
  # expresion list for parameter and jacobian evaluation
  expressionList <- lapply(trafo, function(myrel) parse(text=as.character(myrel)))
  listExpression <- parse(text = paste("list(", paste(trafo, collapse=", "), ")"))
  listJac <- unlist(lapply(parameters, function(var) {
    unlist(lapply(expressionList, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
  }))
  listJac <- parse(text = paste("list(", paste(listJac, collapse=", "), ")"))
  jacNames <- expand.grid(names(trafo), parameters)
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv", exact = TRUE)
    
    # replace fixed parameters by values given in fixed
    if(!is.null(fixed)) {
      is.fixed <- which(names(p)%in%names(fixed)) 
      if(length(is.fixed)>0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    # check for parameters which are not defined in parameters
    emptypars <- names(p)[!names(p)%in%parameters & !names(p)%in%names(fixed)]
    
    # compute transformation output
    out <- with(as.list(p), { 
      
      values <- unlist(eval(listExpression))
      names(values) <- names(trafo)
      
      jacValues <- unlist(eval(listJac))
      jacobian <- matrix(jacValues, ncol=length(parameters), nrow=length(trafo))
      if(any(parameters %in% names(fixed))) {
        jacobian <- matrix(jacobian[,-which(parameters %in% names(fixed))], nrow=length(trafo))
        colnames(jacobian) <- parameters[!parameters%in%names(fixed)]
        rownames(jacobian) <- names(trafo)
      } else {
        colnames(jacobian) <- parameters
        rownames(jacobian) <- names(trafo)
      }
      
      # Append zeros for emptypars
      if(length(emptypars) > 0) {
        empty <- matrix(0, nrow=dim(jacobian)[1], ncol=length(emptypars))
        colnames(empty) <- emptypars
        jacobian <- cbind(jacobian, empty)
      }    
      
      # Multiplication with deriv of p
      if(!is.null(dP)) jacobian <- jacobian%*%dP[colnames(jacobian),]
      
      
      if(deriv) attr(values, "deriv") <- jacobian
      
      return(values)
      
    })
    
    
    return(out)
  }
  
  class(p2p) <- "par"
  return(p2p)
  
}



#' Parameter transformation (implicit)
#' 
#' @param trafo Named character vector defining the equations to be set to zero. 
#' Names correspond to dependent variables.
#' @param parameters Character vector, the independent variables.  
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @details Usually, the equations contain the dependent variables, the independent variables and 
#' other parameters. The argument \code{p} of \code{p2p} must provide values for the independent
#' variables and the parameters but ALSO FOR THE DEPENDENT VARIABLES. Those serve as initial guess
#' for the dependent variables. The dependent variables are then numerically computed by 
#' \link{rootSolve::multiroot}. The Jacobian of the solution with respect to dependent variables
#' and parameters is computed by the implicit function theorem. The function \code{p2p} returns
#' all parameters as they are with corresponding 1-entries in the Jacobian.
#' #' @seealso \link{P} for explicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' ########################################################################
#' ## Example 1: Steady-state trafo
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
#' P.steadyState(p.outerValues)
#' 
#' ########################################################################
#' ## Example 2: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState %o% P.log)(p.outerValue)
Pi <- function(trafo, parameters=NULL) {

  
  
  states <- names(trafo)
  nonstates <- getSymbols(trafo, exclude = states)
  dependent <- setdiff(states, parameters)
  
  trafo.alg <- funC.algebraic(trafo[dependent])
  ftrafo <- function(x, parms) {
    out <- trafo.alg(as.list(c(x, parms)))
    structure(as.numeric(out), names = colnames(out))
  }
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv")
    
    # replace fixed parameters by values given in fixed
    if(!is.null(fixed)) {
      is.fixed <- which(names(p)%in%names(fixed)) 
      if(length(is.fixed)>0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    # check for parameters which are not computed by multiroot
    emptypars <- names(p)[!names(p)%in%c(dependent, fixed)]
    
    # Compute steady state concentrations
    myroot <- rootSolve::multiroot(ftrafo, 
                                   start = p[dependent], 
                                   parms = p[setdiff(names(p), dependent)])
    
    # Output parameters
    out <- c(myroot$root, p[setdiff(names(p), names(myroot$root))])
    
    # Compute jacobian d(root)/dp
    dfdx <- rootSolve::gradient(ftrafo, 
                                x = myroot$root, 
                                parms = p[setdiff(names(p), names(myroot$root))])
    dfdp <- rootSolve::gradient(ftrafo, 
                                x = p[setdiff(names(p), names(myroot$root))],
                                parms = myroot$root)
    dxdp <- solve(dfdx, -dfdp)
    #print(dxdp)
        
       
    # Assemble total jacobian
    jacobian <- matrix(0, length(out), length(p))
    colnames(jacobian) <- names(p)
    rownames(jacobian) <- names(out)
    for(ep in emptypars) jacobian[ep, ep] <- 1
    jacobian[rownames(dxdp), colnames(dxdp)] <- dxdp 
    jacobian <- jacobian[,setdiff(names(p), names(fixed))]
    
    # Multiplication with deriv of p
    if(!is.null(dP)) jacobian <- jacobian%*%dP[colnames(jacobian),]
    
    
    if(deriv) attr(out, "deriv") <- jacobian
    return(out)
    
  }
  
  class(p2p) <- "par"
  return(p2p)
  
}


#' Concatenation of parameter transformations
#' 
#' @param p1 Return value of \link{P} or \link{Pi}
#' @param p2 Return value of \link{P} or \link{Pi}
#' @return A function \code{p2p(p, fixed = NULL, deriv = TRUE)}, the concatenation of \code{p1} and 
#' \code{p2}.
#' @aliases concatenation
#' @examples
#' #' ########################################################################
#' ## Example: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pi(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState %o% P.log)(p.outerValue)
"%o%" <- function(p1, p2) function(p, fixed=NULL, deriv = TRUE) p1(p2(p, fixed, deriv), fixed, deriv)

