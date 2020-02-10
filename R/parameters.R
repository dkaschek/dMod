## Functions to generate parameter transformation ----


#' Generate a parameter transformation function
#' @description  Generate parameter transformation function from a
#' named character vector or object of class \link{eqnvec}. This is a wrapper
#' function for \link{Pexpl} and \link{Pimpl}. See for more details there.
#' @param trafo object of class \code{eqnvec} or named character or list thereof. In case,
#' trafo is a list, \code{P()} is called on each element and conditions are assumed to be
#' the list names.
#' @param parameters character vector
#' @param condition character, the condition for which the transformation is generated
#' @param attach.input attach those incoming parameters to output which are not overwritten by
#' the parameter transformation.
#' @param keep.root logical, applies for \code{method = "implicit"}. The root of the last
#' evaluation of the parameter transformation function is saved as guess for the next 
#' evaluation.
#' @param compile logical, compile the function (see \link{funC0})
#' @param modelname character, see \link{funC0}
#' @param method character, either \code{"explicit"} or \code{"implicit"}
#' @param verbose Print out information during compilation
#' @return An object of class \link{parfn}.
#' @export
P <- function(trafo = NULL, parameters=NULL, condition = NULL, attach.input = FALSE,  keep.root = TRUE, compile = FALSE, modelname = NULL, method = c("explicit", "implicit"), verbose = FALSE) {
  
  if (is.null(trafo)) return()
  if (!is.list(trafo)) {
    trafo <- list(trafo)
    names(trafo) <- condition
  }
  
  method <- match.arg(method)
  
  Reduce("+", lapply(1:length(trafo), function(i) {
  
      switch(method, 
           explicit = Pexpl(trafo = as.eqnvec(trafo[[i]]), parameters = parameters, attach.input = attach.input, condition = names(trafo[i]), compile = compile, modelname = modelname, verbose = verbose),
           implicit = Pimpl(trafo = as.eqnvec(trafo[[i]]), parameters = parameters, keep.root = keep.root, condition = names(trafo[i]), compile = compile, modelname = modelname, verbose = verbose))
    
    
  }))
  
  
}

#' Parameter transformation
#' 
#' @param trafo Named character vector. Names correspond to the parameters being fed into
#' the model (the inner parameters). The elements of tafo are equations that express 
#' the inner parameters in terms of other parameters (the outer parameters)
#' @param parameters Character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(trafo)} the identity transformation is assumed.
#' @param attach.input attach those incoming parameters to output which are not overwritten by
#' the parameter transformation. 
#' @param compile Logical, compile the function (see \link{funC0})
#' @param condition character, the condition for which the transformation is generated
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to R command line.
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @seealso \link{Pimpl} for implicit parameter transformations
#' @examples
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", 
#'               A = "exp(logA)", B = "exp(logB)")
#' p_log <- P(logtrafo)
#' 
#' pars <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' out <- p_log(pars)
#' getDerivs(out)
#' @export
Pexpl <- function(trafo, parameters=NULL, attach.input = FALSE, condition = NULL, compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  # get outer parameters
  symbols <- getSymbols(trafo)
  
  if(is.null(parameters)) {
    parameters <- symbols 
  } else {
    identity <- parameters[which(!parameters%in%symbols)]
    names(identity) <- identity
    trafo <- c(trafo, identity)
  }
  
  # Modify modelname by condition
  if (!is.null(modelname) && !is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  modelname_deriv <- NULL
  if (!is.null(modelname)) modelname_deriv <- paste(modelname, "deriv", sep = "_")

  dtrafo <- jacobianSymb(trafo, parameters)
  
  PEval <- funC0(trafo, parameters = parameters, compile = compile, modelname = modelname, 
                 verbose = verbose, convenient = FALSE, warnings = FALSE)
  dPEval <- funC0(dtrafo, parameters = parameters, compile = compile, 
                  modelname = modelname_deriv, 
                  verbose = verbose, convenient = FALSE, warnings = FALSE)
  
  
  # Controls to be modified from outside
  controls <- list(attach.input = attach.input)
  
  # the parameter transformation function to be returned
  p2p <- function(pars, fixed=NULL, deriv = TRUE) {
    
    p <- pars
    attach.input <- controls$attach.input
    
    # Evaluate transformation
    args <- c(p, fixed)
    pinner <- PEval(p = args)[1,]

    myderiv <- NULL
    if(deriv) {
      
      jac.matrix <- matrix(0, nrow = length(pinner), ncol = length(args), dimnames = list(names(pinner), names(args)))
      jac.matrix[1:length(pinner), match(parameters, names(args))] <- dPEval(p = args)[1,]
      jac.matrix <- jac.matrix[, names(p), drop = FALSE] # delete fixed
  
      dP <- attr(p, "deriv", exact = TRUE)
      if(!is.null(dP)) jac.matrix <- jac.matrix %*% submatrix(dP, rows = colnames(jac.matrix))
      
      myderiv <- jac.matrix
    }
    
    pinner <- as.parvec(pinner, deriv = myderiv)
    if (attach.input & !all(names(pars) %in% names(pinner)))
      pinner <- c(pinner, as.parvec(pars[setdiff(names(pars), names(pinner))]))
    
    return(pinner)
    
  }
  
  attr(p2p, "equations") <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname") <- modelname
  
  parfn(p2p, parameters, condition)
  
}


#' Parameter transformation (implicit)
#' 
#' @param trafo Named character vector defining the equations to be set to zero. 
#' Names correspond to dependent variables.
#' @param parameters Character vector, the independent variables.  
#' @param condition character, the condition for which the transformation is generated
#' @param compile Logical, compile the function (see \link{funC0})
#' @param keep.root logical, applies for \code{method = "implicit"}. The root of the last
#' evaluation of the parameter transformation function is saved as guess for the next 
#' evaluation.
#' @param positive logical, returns projection to the (semi)positive range. Comes with a warning if
#' the steady state has been found to be negative. 
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
#' @seealso \link{Pexpl} for explicit parameter transformations
#' @examples
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
#' (P.steadyState * P.log)(p.outerValue)
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
#' @export
#' @import rootSolve
Pimpl <- function(trafo, parameters=NULL, condition = NULL, keep.root = TRUE, positive = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  
  states <- names(trafo)
  nonstates <- getSymbols(trafo, exclude = states)
  dependent <- setdiff(states, parameters)

  # Modify modelname by condition
  if (!is.null(modelname) && !is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Then add suffix(es) for derivative function
  modelname_dfdx <- NULL
  modelname_dfdp <- NULL
  if (!is.null(modelname)) {
    modelname_dfdx <- paste(modelname, "dfdx", sep = "_")
    modelname_dfdp <- paste(modelname, "dfdp", sep = "_")
  }
  
  
  # Introduce a guess where Newton method starts
  guess <- NULL
  
  trafo0 <- trafo[dependent]
  trafo0.dfdx <- jacobianSymb(trafo0, dependent)
  trafo0.dfdp <- jacobianSymb(trafo0, c(nonstates, parameters))
  
  
  trafo.alg <- funC0(trafo0, parameters = c(states, nonstates), 
                     compile = compile, modelname = modelname, verbose = verbose,
                     convenient = FALSE, warnings = FALSE)
  trafo.alg.dfdx <- funC0(trafo0.dfdx, parameters = c(states, nonstates), 
                          compile = compile, modelname = modelname_dfdx, verbose = verbose,
                          convenient = FALSE, warnings = FALSE)
  trafo.alg.dfdp <- funC0(trafo0.dfdp, parameters = c(states, nonstates), 
                          compile = compile, modelname = modelname_dfdp, verbose = verbose,
                          convenient = FALSE, warnings = FALSE)
  
  ftrafo <- function(x, parms) {
    trafo.alg(p = c(x, parms))[1,]
  }
  ftrafo.dfdx <- function(x, parms) {
    matrix(trafo.alg.dfdx(p = c(x, parms))[1,], nrow = length(dependent), ncol = length(dependent), dimnames = list(dependent, dependent))
  }
  ftrafo.dfdp <- function(x, parms) {
    matrix(trafo.alg.dfdp(p = c(x, parms))[1,], nrow = length(dependent), ncol = length(nonstates) + length(parameters), dimnames = list(dependent, c(nonstates, parameters)))
  }
  
  
  # Controls to be modified from outside
  controls <- list(keep.root = keep.root,
                   positive = positive)
  
  # the parameter transformation function to be returned
  p2p <- function(pars, fixed=NULL, deriv = TRUE) {
    
    
    p <- pars
    keep.root <- controls$keep.root
    positive <- controls$positive
    
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
    
    # Set guess if available
    p0 <- p
    if(!is.null(guess)) 
      p[intersect(dependent, names(guess))] <- guess[intersect(dependent, names(guess))]
    
    # If no initial guesses provides, set to 1
    if (!all(dependent %in% names(p)))
      p[setdiff(dependent, names(p))] <- 1
    
    getRoot <- function(p) {
      
      # Compute steady state concentrations
      rootSolve::multiroot(ftrafo, positive = FALSE,
                           start = p[dependent], 
                           parms = p[setdiff(names(p), dependent)])
    }
    
    myroot <- getRoot(p)
    if (any(myroot$root < 0) & positive) {
      
      p <- p0
      
      # try again with empty guess
      myroot <- getRoot(p)
      
      myroot$root[myroot$root < 0] <- 0
      warning("Found negative steady state. Negative elements have been set to 0.")
      
      out <- c(myroot$root, p[setdiff(names(p), names(myroot$root))])
      if (keep.root) guess <<- NULL
      
    } else {
      
      # Output parameters, write in outer environment if doGuess
      out <- c(myroot$root, p[setdiff(names(p), names(myroot$root))])
      if (keep.root) guess <<- out
      
    }
    
    
    # Compute jacobian d(root)/dp
    dfdx <- ftrafo.dfdx(x = myroot$root, parms = p[setdiff(names(p), names(myroot$root))])
    dfdp <- ftrafo.dfdp(x = p[setdiff(names(p), names(myroot$root))], parms = myroot$root)
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
  attr(p2p, "modelname") <- modelname
  
  parfn(p2p, parameters, condition)
  
  
  
  
  
}




## Functions to simplify the creation of parameter transformations ----

#' Define parameter transformations by \code{define()}, \code{branch()} and \code{insert()}
#' 
#' @param trafo named character vector of parametric expressions or object 
#' of class \code{eqnvec}
#' @param expr character of the form \code{"lhs ~ rhs"} where both \code{lhs}
#' and \code{rhs} can contain a number of symbols for which vaues are passed
#' by the \code{...} argument
#' @param  conditionMatch optional character, Use as regular expression to apply the reparameterization only to conditions containing conditionMatch
#' @param ... used to pass values for symbols as named arguments
#' @return object of the same class as trafo or list thereof, if \code{branch()} has been
#' used.
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(mytrafo[[i]])
    
    with(mytable, {
      args <- c(list(expr = expr, trafo = mytrafo[[i]], reset = TRUE), eval(dots))
      do.call(repar, args)
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(.currentTrafo)
    
    
    
    with(mytable, {
      .fun <- function() {
        # subset conditions by logicals expressions supplied by the dots
        dots_eval <- eval(dots)                                            # convert from substituted to language
        if (length(dots_eval) == 0) 
          return(do.call(repar, list(expr = expr, trafo = .currentTrafo)))
        
        dots_eval_eval <- lapply(dots_eval, function(i) eval.parent(i, 3)) # evaluate the language in the "mytable" frame. parent1: lapply, parent2: .fun, parent3: with
        which_logical <- vapply(dots_eval_eval, function(i) {is.logical(i)}, FUN.VALUE = vector("logical", 1)) # which of the dots are logical
        logical_dots <- dots_eval[which_logical]                           # subset to matching conditions
        matching <- do.call(c, logical_dots)
        if(!is.null(matching)) { # null means no logical dots were supplied
          if (any(!matching)) {  
            return(.currentTrafo) }}
        
        args <- c(list(expr = expr, trafo = .currentTrafo), dots_eval_eval[!which_logical]) # feed the rest of the eval'd dots into repar
        return(do.call(repar, args))
      }
      .fun()
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}





#' @export
#' @rdname define
#' @param table table of covariates as data frame. Rownames are used as unique identifier,
#' usually called "conditions", and columns represent covariates associated with these conditions.
#' @param conditions character vector with condition names. Overwrites the rownames of table.
branch <- function(trafo, table = NULL, conditions = rownames(table)) {
  
  if (is.null(table) & is.null(conditions)) 
    return(trafo)
  
  if (is.null(conditions)) conditions <- paste0("C", 1:nrow(table))
  if (is.null(table)) table <- data.frame(condition = conditions, row.names = conditions)
  
  out <- lapply(conditions, function(x) trafo)
  names(out) <- conditions
  
  rownames(table) <- conditions
  attr(out, "tree") <- table
  
  return(out)  
  
}


#' Reparameterization
#' 
#' @param expr character of the form \code{"lhs ~ rhs"} where \code{rhs}
#' reparameterizes \code{lhs}. Both \code{lhs} and \code{rhs}
#' can contain a number of symbols whose values need to be passed by the \code{...} argument.
#' @param trafo character or equation vector or list thereof. The object where the replacement takes place in
#' @param ... pass symbols as named arguments
#' @param reset logical. If true, the trafo element corresponding to lhs is reset according to rhs. 
#' If false, lhs wherever it occurs in the rhs of trafo is replaced by rhs of the formula.
#' @return an equation vector with the reparameterization.
#' @details Left and right-hand side of \code{expr} are searched for symbols. If separated by
#' "_", symbols are recognized as such, e.g. in \code{Delta_x} where the symbols are 
#' "Delta" and "x". Each symbol for which values (character or numbers) are passed by the
#' \code{...} argument is replaced.
#' @export
#' @importFrom stats as.formula
#' @examples
#' innerpars <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' 
#' trafo <- repar("x ~ x", x = innerpars)
#' trafo <- repar("x ~ y", trafo, x = names(constraints), y = constraints)
#' trafo <- repar("x ~ exp(x)", trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition", trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ..., reset = FALSE) {
  
  if (inherits(expr, "formula")) expr <- deparse(expr)
  
  parsed.expr <- as.character(stats::as.formula(gsub("_", ":", expr, fixed = TRUE)))
  lhs <- parsed.expr[2]
  lhs.symbols <- getSymbols(lhs)
  rhs <- parsed.expr[3]
  rhs.symbols <- getSymbols(rhs)
  
  # Make sure that arguments are characters
  args <- lapply(list(...), as.character)
  
  replacements <- as.data.frame(args, stringsAsFactors = FALSE)
  
  lhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], lhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  rhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], rhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  if (is.null(trafo)) {
    trafo <- as.eqnvec(structure(lhs, names = lhs))
  } else if (is.list(trafo) & !reset) {
    trafo <- lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  } else if (is.character(trafo) & !reset) {
    trafo <- replaceSymbols(lhs, rhs, trafo)
  } else if (is.list(trafo) & reset) {
    trafo <- lapply(trafo, function(t) {t[lhs] <- rhs; return(t)})
  } else if (is.character(trafo) & reset) {
    trafo[lhs] <- rhs
  }
  
  return(trafo)
  
  
}

paste_ <- function(...) paste(..., sep = "_")



