#' Generate a parameter transformation function
#' @description  Generate parameter transformation function from a
#' named character vector or object of class \link{eqnvec}. This is a wrapper
#' function for \link{Pexpl} and \link{Pimpl}. See for more details there.
#' @param trafo object of class \code{eqnvec} or named character
#' @param parameters character vector
#' @param compile logical
#' @param modelname character
#' @param method character, either \code{"explicit"} or \code{"implicit"}
#' @param verbose Print out information during compilation
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)}
#' @export
P <- function(trafo = NULL, parameters=NULL, condition = NULL, compile = FALSE, modelname = NULL, method = c("explicit", "implicit"), verbose = FALSE) {
  
  if (is.null(trafo)) return()
  
  method <- match.arg(method)
  switch(method, 
         explicit = Pexpl(trafo = trafo, parameters = parameters, condition = condition, compile = compile, modelname = modelname, verbose = verbose),
         implicit = Pimpl(trafo = trafo, parameters = parameters, condition = condition, compile = compile, modelname = modelname, verbose = verbose))
  
}

#' The identity parameter transformation
#' @export
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} returning \code{p}
#' with unit matrix as derivative or derivative inherited from \code{p}.
P0 <- function() {
  
  myfn <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv", exact = TRUE)
    
    
    myderiv <- NULL
    if (deriv) {
      myderiv <- diag(x = 1, nrow = length(p), ncol = length(p))
      rownames(myderiv) <- colnames(myderiv) <- names(p)
      if (!is.null(dP)) myderiv <- myderiv %*% submatrix(dP, rows = colnames(myderiv))
    }
    
    as.parvec(p, deriv = myderiv)
    
  }
  
  return(myfn)
    
}

#' Parameter transformation
#' 
#' @param trafo Named character vector. Names correspond to the parameters being fed into
#' the model (the inner parameters). The elements of tafo are equations that express 
#' the inner parameters in terms of other parameters (the outer parameters)
#' @param parameters Character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(trafo)} the identity transformation is assumed.
#' @param compile Logical, compile the function (see \link{funC0})
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to R command line.
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @seealso \link{Pi} for implicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' \dontrun{
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' }
#' @export
Pexpl <- function(trafo, parameters=NULL, condition = NULL, compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  # get outer parameters
  symbols <- getSymbols(trafo)
  
  if(is.null(parameters)) {
    parameters <- symbols 
  } else {
    identity <- parameters[which(!parameters%in%symbols)]
    names(identity) <- identity
    trafo <- c(trafo, identity)
  }
  
  
  # expression list for parameter and jacobian evaluation
  #trafo.list <- lapply(trafo, function(myrel) parse(text=as.character(myrel)))
  #jacobian <- unlist(lapply(parameters, function(var) {
  #  unlist(lapply(trafo.list, function(myexp) paste(deparse(D(myexp, as.character(var))), collapse="")))
  #}))
  dtrafo <- jacobianSymb(trafo, parameters)
  
  #jacNames <- expand.grid.alt(names(trafo), parameters)
  #jacNames <- paste(jacNames[,1], jacNames[,2], sep=".")
  
  #dtrafo <- jacobian; names(dtrafo) <- jacNames
  
  PEval <- funC0(trafo, compile = compile, modelname = modelname, verbose = verbose)
  dPEval <- funC0(dtrafo, compile = compile, modelname = paste(modelname, "deriv", sep = "_"), verbose = verbose)
  
  # the parameter transformation function to be returned
  p2p <- function(p, fixed=NULL, deriv = TRUE) {
    
    # Inherit from p
    dP <- attr(p, "deriv", exact = TRUE)
    
    # Evaluate transformation
    args <- c(as.list(p), as.list(fixed))
    pinner <- PEval(args)[1,]
    dpinner <- dPEval(args)[1,]
    
    # Construct output jacobian
    jac.vector <- rep(0, length(pinner)*length(p))
    names(jac.vector) <- outer(names(pinner), names(p), function(x, y) paste(x, y, sep = "."))
    
    names.intersect <- intersect(names(dpinner), names(jac.vector))
    jac.vector[names.intersect] <- as.numeric(dpinner[names.intersect])
    jac.matrix <- matrix(jac.vector, length(pinner), length(p), dimnames = list(names(pinner), names(p)))
    
    
    if(!is.null(dP)) jac.matrix <- jac.matrix%*%submatrix(dP, rows = colnames(jac.matrix))
    
    myderiv <- NULL
    if(deriv) myderiv <- jac.matrix
    
    as.parvec(pinner, deriv = myderiv)
    
  }
  
  attr(p2p, "equations") <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  
  parfn(p2p, parameters, condition)
  
}


#' Parameter transformation (implicit)
#' 
#' @param trafo Named character vector defining the equations to be set to zero. 
#' Names correspond to dependent variables.
#' @param parameters Character vector, the independent variables.  
#' @param compile Logical, compile the function (see \link{funC0})
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to R command line.
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @details Usually, the equations contain the dependent variables, the independent variables and 
#' other parameters. The argument \code{p} of \code{p2p} must provide values for the independent
#' variables and the parameters but ALSO FOR THE DEPENDENT VARIABLES. Those serve as initial guess
#' for the dependent variables. The dependent variables are then numerically computed by 
#' \link[rootSolve]{multiroot}. The Jacobian of the solution with respect to dependent variables
#' and parameters is computed by the implicit function theorem. The function \code{p2p} returns
#' all parameters as they are with corresponding 1-entries in the Jacobian.
#' @seealso \link{Pexpl} for explicit parameter transformations and
#' \link{concatenation} for the concatenation of parameter transformations
#' @examples
#' \dontrun{
#' ########################################################################
#' ## Example 1: Steady-state trafo
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#' 
#' p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
#' P.steadyState(p.outerValues)
#' 
#' ########################################################################
#' ## Example 2: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState %o% P.log)(p.outerValue)
#' 
#' ########################################################################
#' ## Example 3: Steady-states with conserved quantitites
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B", B = "k1*A - k2*B")
#' replacement <- c(B = "A + B - total")
#' f[names(replacement)] <- replacement
#' 
#' pSS <- Pimpl(f, "total")
#' pSS(c(k1 = 1, k2 = 2, A = 5, B = 5, total = 3))
#' }
#' @export
Pimpl <- function(trafo, parameters=NULL, condition = NULL, keep.root = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE) {

  states <- names(trafo)
  nonstates <- getSymbols(trafo, exclude = states)
  dependent <- setdiff(states, parameters)
  
  # Introduce a guess where Newton method starts
  guess <- NULL
  
  trafo.alg <- funC0(trafo[dependent], compile = compile, modelname = modelname, verbose = verbose)
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

    # Set guess
    if(!is.null(guess)) 
      p[intersect(dependent, names(guess))] <- guess[intersect(dependent, names(guess))]
    
    # check for parameters which are not computed by multiroot
    emptypars <- names(p)[!names(p)%in%c(dependent, fixed)]
    
    # Compute steady state concentrations
    myroot <- rootSolve::multiroot(ftrafo, positive = TRUE,
                                   start = p[dependent], 
                                   parms = p[setdiff(names(p), dependent)])
    
    # Output parameters, write in outer environment if doGuess
    out <- c(myroot$root, p[setdiff(names(p), names(myroot$root))])
    if(keep.root) guess <<- out
    
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
    jacobian <- submatrix(jacobian, cols = setdiff(names(p), names(fixed)))
    
    # Multiplication with deriv of p
    if(!is.null(dP)) jacobian <- jacobian%*%submatrix(dP, rows = colnames(jacobian))
    
    myderiv <- NULL
    if(deriv) myderiv <- jacobian
    
    as.parvec(out, deriv = myderiv)
    
  }
  
  attr(p2p, "equations") <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
    
  parfn(p2p, parameters, condition)
  

  

  
}

#' Concatenation of parameter transformations
#' 
#' @param p1 Return value of \link{P} or \link{Pi}
#' @param p2 Return value of \link{P} or \link{Pi}
#' @return A function \code{p2p(p, fixed = NULL, deriv = TRUE)}, the concatenation of \code{p1} and 
#' \code{p2}.
#' @aliases concatenation
#' @examples
#' \dontrun{
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
#' }
#' @export
"%o%" <- function(p1, p2) function(p, fixed=NULL, deriv = TRUE) p1(p2(p, fixed = fixed, deriv = deriv), deriv = deriv)




