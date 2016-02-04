## ODE model class -------------------------------------------------------------------


#' Generate the model objects for use in Xs (models with sensitivities)
#' 
#' @param f Named character vector with the ODE
#' @param deriv logical, generate sensitivities or not
#' @param forcings Character vector with the names of the forcings
#' @param fixed Character vector with the names of parameters (initial values and dynamic) for which
#' no sensitivities are required (will speed up the integration).
#' @param modelname Character, the name of the C file being generated.
#' @param verbose Print compiler output to R command line.
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object)
#' @export
#' @import cOde
odemodel <- function(f, deriv = TRUE, forcings=NULL, fixed=NULL, modelname = "odemodel", verbose = FALSE, ...) {
  
  modelname_s <- paste0(modelname, "_s")
  
  func <- cOde::funC(f, forcings = forcings, modelname = modelname , ...)
  extended <- NULL
  if (deriv) {  
    s <- sensitivitiesSymb(f, 
                           states = setdiff(attr(func, "variables"), fixed), 
                           parameters = setdiff(attr(func, "parameters"), fixed), 
                           inputs = forcings,
                           reduce = TRUE)
    fs <- c(f, s)
    outputs <- attr(s, "outputs")
    extended <- cOde::funC(fs, forcings = forcings, outputs = outputs, modelname = modelname_s, ...)
  }  
  
  out <- list(func = func, extended = extended)
  attr(out, "class") <- "odemodel"
  return(out)
  
  
}


## Equation classes -------------------------------------------------------

#' Generate equation vector object
#'
#' @description The eqnvec object stores explicit algebraic equations, like the
#' right-hand sides of an ODE, observation functions or parameter transformations
#' as named character vectors.
#' @param ... mathematical expressions as characters to be coerced,
#' the right-hand sides of the equations
#' @return object of class \code{eqnvec}, basically a named character.
#' @export
eqnvec <- function(...) {
  
  mylist <- list(...)
  if (length(mylist) > 0) {
    mynames <- paste0("eqn", 1:length(mylist))
    is.available <- !is.null(names(mylist))
    mynames[is.available] <- names(mylist)[is.available]  
    
    names(mylist) <- mynames
    out <- unlist(mylist)
    
    return(as.eqnvec(out))
    
  } else {
    
    return(NULL)
  
  }
    
}

#' Generate eqnlist object
#'
#' @description The eqnlist object stores an ODE as a list of stoichiometric matrix,
#' rate expressions, state names and compartment volumes.
#' @export
#' @param smatrix Matrix of class numeric. The stoichiometric matrix,
#' one row per reaction/process and one column per state.
#' @param states Character vector. Names of the states.
#' @param rates Character vector. The rate expressions.
#' @param volumes Named character, volume parameters for states. Names must be a subset of the states.
#' Values can be either characters, e.g. "V1", or numeric values for the volume. If \code{volumes} is not
#' \code{NULL}, missing entries are treated as 1.
#' @param description Character vector. Description of the single processes.
#' @return An object of class \code{eqnlist}, basically a list.
eqnlist <- function(smatrix = NULL, states = colnames(smatrix), rates = NULL, volumes = NULL, description = NULL) {

  # Dimension checks and preparations for non-empty argument list.
  if (all(!is.null(c(smatrix, states, rates)))) {
    #Dimension checks
    d1 <- dim(smatrix)
    l2 <- length(states)
    l3 <- length(rates)
    if (l2 != d1[2]) stop("Number of states does not coincide with number of columns of stoichiometric matrix")
    if (l3 != d1[1]) stop("Number of rates does not coincide with number of rows of stoichiometric matrix")

    # Prepare variables
    smatrix <- as.matrix(smatrix)
    colnames(smatrix) <- states
    if (is.null(description)) {
      description <- 1:nrow(smatrix)
    }
  }

  out <- list(smatrix = smatrix,
              states = as.character(states),
              rates = as.character(rates),
              volumes = volumes,
              description = as.character(description))
  class(out) <- c("eqnlist", "list")

  return(out)
}



## Parameter classes --------------------------------------------------------

#' @export
parfn <- function(p2p, parameters = NULL, condition = NULL) {
  
  force(condition)
  mappings <- list()
  mappings[[1]] <- p2p
  names(mappings) <- condition
  
  outfn <- function(p, fixed = NULL, deriv = TRUE, conditions = condition) {
    
    
    result <- p2p(p = p, fixed = fixed, deriv = deriv, ...)
    # If NULL condition, valid for all conditions
    # else feed into correct slot 
    if (is.null(condition)) {
      outlist <- lapply(conditions, function(i) result)
      names(outlist) <- conditions
    } else {
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      outlist[[condition]] <- result
    }
    
    return(outlist)
    
  }
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- condition
  class(outfn) <- c("parfn", "fn")
  return(outfn)
  
  
}

#' Generate a paramter frame
#'
#' @description A parameter frame is a data.frame where the rows correspond to different
#' parameter specifications. The columns are divided into three parts. (1) the meta-information
#' columns (e.g. index, value, constraint, etc.), (2) the attributes of an objective function
#' (e.g. data contribution and prior contribution) and (3) the parameters.
#' @param x data.frame.
#' @param parameters character vector, the names of the parameter columns.
#' @param metanames character vector, the names of the meta-information columns.
#' @param obj.attributes character vector, the names of the objective function attributes.
#' @return An object of class \code{parframe}, i.e. a data.frame with attributes for the
#' different names. Inherits from data.frame.
#' @export
parframe <- function(x = NULL, parameters = colnames(x), metanames = NULL, obj.attributes = NULL) {

  if (!is.null(x)) {
    rownames(x) <- NULL
    out <- as.data.frame(x)
  } else {
    out <- data.frame()
  }

  attr(out, "parameters") <- parameters
  attr(out, "metanames") <- metanames
  attr(out, "obj.attributes") <- obj.attributes
  class(out) <- c("parframe", "data.frame")

  return(out)

}

#' Parameter list
#' 
#' @description The special use of a parameter list is to save
#' the outcome of multiple optimization runs, e.g. by \code{mstrust},
#' into one list.#' 
#' @param ... Objects to be coerced to parameter list.
#' @export
parlist <- function(...) {
  
  mylist <- list(...)
  return(as.parlist(mylist))
  
}



#' Parameter vector
#'
#' @description A parameter vector is a named numeric vector (the parameter values)
#' together with a "deriv" attribute (the Jacobian of a parameter transformation by which
#' the parameter vector was generated).
#' @param deriv matrix with rownames (according to names of \code{...}) and colnames
#' according to the names of the parameter by which the parameter vector was generated.
#' @return An object of class \code{parvec}, i.e. a named numeric vector with attribute "deriv".
#'
#' @export
parvec <- function(..., deriv = NULL) {

  mylist <- list(...)
  if (length(mylist) > 0) {
    mynames <- paste0("par", 1:length(mylist))
    is.available <- !is.null(names(mylist))
    mynames[is.available] <- names(mylist)[is.available]  
    
    out <- as.numeric(unlist(mylist))
    names(out) <- mynames
    
    return(as.parvec(out, deriv = deriv))
    
  } else {
    
    return(NULL)
    
  }
  

}


## Prediction classes ----------------------------------------------------

#' Prediction function
#'
#' @description A prediction function is a function \code{x(times, pars, fixed, deriv, ...)}
#' which returns a list of model predictions. Each entry of the list is a \link{prdframe}
#' as being produced by a low-level prediction function, see e.g. \link{Xs}. The different entries
#' correspond to different experimental conditions which are supposed to be matched to a
#' \link{datalist}.
#' @return Object of class \code{prdfn}, i.e. a function \code{x(times, pars, fixed, deriv, ...)}
#' which returns a \link{prdlist}.
#' @export
prdfn <- function(P2X, parameters = NULL, condition = NULL) {
  
  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- P2X
  names(mappings) <- condition
  
  outfn <- function(times, pars, fixed = NULL, deriv = TRUE, conditions = mycondition) {
    
    
    result <- P2X(times = times, pars = pars, deriv = deriv)
    # If NULL condition, valid for all conditions
    # else feed into correct slot 
    if (is.null(condition)) {
      outlist <- lapply(conditions, function(i) result)
      names(outlist) <- conditions
    } else {
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      outlist[[condition]] <- result
    }
    
    outlist <- as.prdlist(outlist)
    attr(outlist, "pars") <- pars
    
    return(outlist)
    
  }
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- mycondition
  class(outfn) <- c("prdfn", "fn") 
  return(outfn)
  
}

#' @export
obsfn <- function(X2Y, parameters = NULL, condition = NULL) {
  
  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- X2Y
  names(mappings) <- condition
  
  outfn <- function(out, pars, fixed = NULL, deriv = TRUE, conditions = mycondition) {
    
    result <- X2Y(out = out, pars = pars)
    # If NULL condition, valid for all conditions
    # else feed into correct slot 
    if (is.null(condition)) {
      outlist <- lapply(conditions, function(i) result)
      names(outlist) <- conditions
    } else {
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      outlist[[condition]] <- result
    }
    
    return(outlist)
    
  }
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- mycondition
  class(outfn) <- c("obsfn", "fn")
  return(outfn)
  
}


#' Prediction frame
#'
#' @description A prediction frame is used to store a model prediction in a matrix. The columns
#' of the matrix are "time" and one column per state. The prediction frame has attributes "deriv",
#' the matrix of sensitivities with respect to "outer parameters" (see \link{P}), an attribute
#' "sensitivities", the matrix of sensitivities with respect to the "inner parameters" (the model
#' parameters, left-hand-side of the parameter transformation) and an attributes "parameters", the
#' names of the outer parameters.
#' @param prediction matrix of model prediction
#' @param deriv matrix of sensitivities wrt outer parameters
#' @param sensitivities matrix of sensitivitie wrt inner parameters
#' @param parameters names of the outer paramters
#' @return Object of class \code{prdframe}, i.e. a matrix with other matrices and vectors as attributes.
#' @export
prdframe <- function(prediction = NULL, deriv = NULL, sensitivities = NULL, parameters = NULL) {

  out <- if (!is.null(prediction)) as.matrix(prediction) else matrix()

  attr(out, "deriv") <- deriv
  attr(out, "sensitivities") <- sensitivities
  attr(out, "parameters") <- parameters
  class(out) <- c("prdframe", "matrix")

  return(out)

}

#' Prediction list
#'
#' @description A prediction list is used to store a list of model predictions
#' from different prediction functions or the same prediction function with different
#' parameter specifications. Each entry of the list is a \link{prdframe}.
#' @param ... objects of class \link{prdframe}
#' conditions.
#' @export
prdlist <- function(...) {
  mylist <- list(...)
  mynames <- names(mylist)
  if (is.null(mynames)) mynames <- as.character(1:length(mylist))
  as.prdlist(mylist, mynames)
}



## Data classes ----------------------------------------------------------------

#' Generate a datalist object
#'
#' @description The datalist object stores time-course data in a list of data.frames.
#' The names of the list serve as identifiers, e.g. of an experimental condition, etc.
#' @param ... data.frame objects to be coerced into a list
#' @param dataframe data.frame with additional columns by which a splitting into 
#' a list of data.frames is performed.
#' @param split.by character vector, interaction of these columns is defined as
#' new column by which the split-up is performed.
#' @param mylist list of data.frame, each data.frame is expected to have columns
#' "name" (factor or character),
#' "time" (numeric), "value" (numeric) and "sigma" (numeric).
#' @param mynames character vector of the length of mylist.
#' @return Object of class \code{datalist}.
#' @export
datalist <- function(...) {
  mylist <- list(...)
  mynames <- names(mylist)
  if (is.null(mynames)) mynames <- as.character(1:length(mylist))
  as.datalist(mylist, mynames)
}


## Objective classes ---------------------------------------------------------

#' Objective function
#'
#' @description An objective function is a function \code{obj(pouter, fixed, deriv, ...)} which
#' returns an \link{objlist}. This function is supposed to be used in an optimizer like \code{trust}
#' from the \code{trust} package.
#' @param ... R code expressing objectives additional to the weighted residual sum of squares of
#' data and model prediction. Available keywords are
#' \code{pouter}, \code{fixed} and \code{deriv}. Any other object being used in the
#' expression can either be passed by the \code{...} argument of the returned function or
#' must be available in the global environment.
#' @param data object of class \link{datalist}
#' @param x object of class \link{prdfn}
#' @param pouter named numeric, pre-initialized parameter vector
#' @param conditions character vector, names of the conditions to be evaluated
#' @return Object of class \code{objfn}, i.e. a function \code{obj(pouter, fixed, deriv, ...)}
#' which returns an \link{objlist}.
#' @export
objfn <- function(data, x) {

  mydata <- data
  myx <- x
  timesD <- sort(unique(c(0, do.call(c, lapply(mydata, function(d) d$time)))))

  myfn <- function(pouter, fixed = NULL, deriv=TRUE, conditions = names(data)) {
    
    data <- mydata[conditions]
    x <- myx
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)

    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) wrss(res(mydata[[cn]], prediction[[cn]])))
    out.data <- Reduce("+", out.data)

    # Combine contributions and attach attributes
    out <- out.data
    attr(out, "data") <- out.data$value
    
    return(out)


  }
  class(myfn) <- c("objfn", "fn")
  return(myfn)

}


#' Generate objective list
#'
#' @description An objective list contains an objective value, a gradient, and a Hessian matrix
#' @param value numeric of length 1
#' @param gradient named numeric
#' @param hessian matrix with rownames and colnames according to gradient names
#' @return Object of class \code{objlist}
#' @export
objlist <- function(value, gradient, hessian) {

  out <- list(value = value, gradient = gradient, hessian = hessian)
  class(out) <- c("objlist", "list")
  return(out)

}

#' Objective frame
#'
#' @description An objective frame is supposed to store the residuals of a model prediction
#' with respect to a data frame.
#' @param mydata data.frame as being generated by \link{res}.
#' @param deriv matrix of the derivatives of the residuals with respect to parameters.
#' @return An object of class \code{objframe}, i.e. a data frame with attribute "deriv".
#' @export
objframe <- function(mydata, deriv = NULL) {

  # Check column names
  mydata <- as.data.frame(mydata)
  correct.names <- c("time", "name", "value", "prediction",
                     "sigma", "residual", "weighted.residual")

  ok <- all(correct.names %in% names(mydata))
  if (!ok) stop("mydata does not have required names")

  out <- mydata[, correct.names]
  attr(out, "deriv") <- deriv
  class(out) <- c("objframe", "data.frame")

  return(out)


}




## General concatenation of functions ------------------------------------------

#' @export
"+.fn" <- function(x1, x2) {
  
  if (is.null(x1)) return(x2)
   
  # prdfn + prdfn
  if (inherits(x1, "prdfn") & inherits(x2, "prdfn")) {
    
    mappings.x1 <- attr(x1, "mappings")
    mappings.x2 <- attr(x2, "mappings")
    
    if (is.null(names(mappings.x1)) || is.null(names(mappings.x2))) stop("General functions (NULL names) cannot be coerced.")
    
    
    mappings <- c(mappings.x1, mappings.x2)
    
    outfn <- function(times, pars, fixed = NULL, deriv = TRUE, conditions = names(mappings)) {
      
      available <- intersect(names(mappings), conditions)
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      outpars <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](times = times, pars = pars, deriv = deriv)
        outpars[[C]] <- attr(outlist[[C]], "pars")
        attr(outlist[[C]], "pars") <- NULL
      }
      
      out <- as.prdlist(outlist)
      attr(out, "pars") <- outpars
      return(out)
      
    }
    
    attr(outfn, "mappings") <- mappings
    attr(outfn, "parameters") <- union(attr(x1, "parameters"), attr(x2, "parameters"))
    attr(outfn, "conditions") <- c(attr(x1, "conditions"), attr(x2, "conditions"))
    class(outfn) <- c("prdfn", "fn")
    
    return(outfn)
    
    
  }
  
  # parfn + parfn
  if (inherits(x1, "parfn") & inherits(x2, "parfn")) {
    
    mappings.x1 <- attr(x1, "mappings")
    mappings.x2 <- attr(x2, "mappings")
    
    if (is.null(names(mappings.x1)) || is.null(names(mappings.x2))) stop("General transformations (NULL names) cannot be coerced.")
    
    
    mappings <- c(mappings.x1, mappings.x2)
    
    outfn <- function(p, fixed = NULL, deriv = TRUE, conditions = names(mappings)) {
      
      available <- intersect(names(mappings), conditions)
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](p = p, fixed = fixed, deriv = deriv)
      }
      
      return(outlist)
      
    }
    
    attr(outfn, "mappings") <- mappings
    attr(outfn, "parameters") <- union(attr(x1, "parameters"), attr(x2, "parameters"))
    attr(outfn, "conditions") <- c(attr(x1, "conditions"), attr(x2, "conditions"))
    class(outfn) <- c("parfn", "fn")
    
    return(outfn)
    
  }
  
}

#' @export
"*.fn" <- function(p1, p2) {

  # obsfn * obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "obsfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    
    outfn <- function(out, pars, fixed = NULL, deriv = TRUE, conditions = conditions.p2) {
      
      step1 <- p2(out = out, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      pars1 <- attr(step1, "pars")
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = pars1[[i]], fixed = fixed, deriv = deriv, condition = names(step1)[i])))
      
      out <- as.prdlist(step2)
      attr(out, "pars") <- pars1
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- attr(p2, "conditions")
    class(outfn) <- c("obsfn", "fn")
    
    return(outfn)
    
  }  
  
  
  # obsfn * prdfn
  if (inherits(p1, "obsfn") & inherits(p2, "prdfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    
    outfn <- function(times, pars, fixed = NULL, deriv = TRUE, conditions = conditions.p2) {
      step1 <- p2(times = times, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      pars1 <- attr(step1, "pars")
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = pars1[[i]], fixed = fixed, deriv = deriv, condition = names(step1)[i])))
      
      out <- as.prdlist(step2)
      attr(out, "pars") <- pars1
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- attr(p2, "conditions")
    class(outfn) <- c("prdfn", "fn")
    
    return(outfn)
    
  }  
  
  
  # prdfn * parfn
  if (inherits(p1, "prdfn") & inherits(p2, "parfn")) {
    
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    
    outfn <- function(times, pars, fixed = NULL, deriv = TRUE, conditions = conditions.p2) {
      
      step1 <- p2(p = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(times = times, pars = step1[[i]], deriv = deriv, condition = names(step1)[i])))
      
      out <- as.prdlist(step2)
      attr(out, "pars") <- step1
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- attr(p2, "conditions")
    class(outfn) <- c("prdfn", "fn")
    
    return(outfn)
    
  }
  
  
  # parfn * parfn
  if (inherits(p1, "parfn") & inherits(p2, "parfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    
    outfn <- function(p, fixed=NULL, deriv = TRUE, conditions = conditions.p2) {
      
      step1 <- p2(p = p, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(p = step1[[i]], deriv = deriv, condition = names(step1)[i])))
      return(step2)
      
    }
    
    attr(outfn, "mappings") <- attr(p2, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- attr(p2, "conditions")
    class(outfn) <- c("parfn", "fn")
    
    return(outfn)
    
  }
  


}


