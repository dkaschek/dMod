## ODE model class -------------------------------------------------------------------


#' Generate the model objects for use in Xs (models with sensitivities)
#' 
#' @param f Something that can be converted to \link{eqnvec}, 
#' e.g. a named character vector with the ODE
#' @param deriv logical, generate sensitivities or not
#' @param forcings Character vector with the names of the forcings
#' @param fixed Character vector with the names of parameters (initial values and dynamic) for which
#' no sensitivities are required (will speed up the integration).
#' @param modelname Character, the name of the C file being generated.
#' @param verbose Print compiler output to R command line.
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object)
#' @export
#' @example inst/examples/odemodel.R
#' @import cOde
odemodel <- function(f, deriv = TRUE, forcings=NULL, fixed=NULL, modelname = "odemodel", verbose = FALSE, ...) {
  
  f <- as.eqnvec(f)
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

## Function classes ------------------------------------------------------

match.fnargs <- function(arglist, choices) {
  
  # Catch the case of names == NULL
  if (is.null(names(arglist))) names(arglist) <- rep("", length(arglist))
  
  # exlude named arguments which are not in choices
  arglist <- arglist[names(arglist) %in% c(choices, "")]
  
  # determine available arguments
  available <- choices %in% names(arglist)
  
  if (!all(available)) names(arglist)[names(arglist) == ""] <- choices[!available]
  
  if (any(duplicated(names(arglist)))) stop("duplicate arguments in prdfn/obsfn/parfn function call")
  
  mapping <- match(choices, names(arglist))
  return(mapping)
  
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
#' @example inst/examples/eqnvec.R
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

#' Parameter transformation function
#' 
#' Generate functions that transform one parameter vector into another
#' by means of a transformation, pushing forward the jacobian matrix
#' of the original parameter. 
#' Usually, this function is called internally, e.g. by \link{P}.
#' However, you can use it to add your own specialized parameter
#' transformations to the general framework.
#' @param p2p a transformation function for one condition, i.e. a function
#' \code{p2p(p, fixed, deriv)} which translates a parameter vector \code{p}
#' and a vector of fixed parameter values \code{fixed} into a new parameter
#' vector. If \code{deriv = TRUE}, the function should return an attribute
#' \code{deriv} with the Jacobian matrix of the parameter transformation.
#' @param parameters character vector, the parameters accepted by the function
#' @param condition character, the condition for which the transformation is defined
#' @return object of class \code{parfn}. Contains attributes "mappings", a list of \code{p2p}
#' functions, "parameters", the union of parameters acceted by the mappings and
#' "conditions", the total set of conditions.
#' @seealso \link{sumfn}, \link{P}
#' @export
parfn <- function(p2p, parameters = NULL, condition = NULL) {
  
  force(condition)
  mappings <- list()
  mappings[[1]] <- p2p
  names(mappings) <- condition
  
  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition) {
   
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pars <- arglist[[1]]
    
    
    overlap <- intersect(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap
    
    if (is.null(overlap) | length(overlap) > 0)
      result <- p2p(pars = pars, fixed = fixed, deriv = deriv)
    else
      result <- NULL
    
    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- lapply(1:length.out, function(i) result)
    names(outlist) <- conditions

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
#' @seealso \link{profile}, \link{mstrust}
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
#' the outcome of multiple optimization runs provided by \link{mstrust},
#' into one list. 
#' @param ... Objects to be coerced to parameter list.
#' @export
#' @seealso \link{load.parlist}, \link{plot.parlist}
parlist <- function(...) {
  
  mylist <- list(...)
  return(as.parlist(mylist))
  
}



#' Parameter vector
#'
#' @description A parameter vector is a named numeric vector (the parameter values)
#' together with a "deriv" attribute (the Jacobian of a parameter transformation by which
#' the parameter vector was generated).
#' @param ... objects to be concatenated
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
#' @description A prediction function is a function \code{x(times, pars, fixed, deriv, conditions)}
#' which returns a list of model predictions. Each entry of the list is a \link{prdframe}
#' as being produced by a low-level prediction function, see e.g. \link{Xs}. The different entries
#' correspond to different experimental conditions which can be compared to a
#' \link{datalist}.
#' @param P2X transformation function as being produced by \link{Xs}.
#' @param parameters character vector with parameter names
#' @param condition character, the condition name
#' @details Prediction functions can be "added" by the "+" operator, see \link{sumfn}. Thereby,
#' predictions for different conditions are merged or, overwritten. Prediction functions can
#' also be concatenated with other functions, e.g. observation functions (\link{obsfn}) or 
#' parameter transformation functions (\link{parfn}) by the "*" operator, see \link{prodfn}.
#' @return Object of class \code{prdfn}, i.e. a function \code{x(times, pars, fixed, deriv, conditions)}
#' which returns a \link{prdlist}.
#' @export
prdfn <- function(P2X, parameters = NULL, condition = NULL) {
  
  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- P2X
  names(mappings) <- condition
  
  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = mycondition) {
   
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
    times <- arglist[[1]]
    pars <- arglist[[2]]
    
    
    overlap <- intersect(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap
    
    if (is.null(overlap) | length(overlap) > 0)
      result <- P2X(times = times, pars = pars, deriv = deriv)
    else
      result <- NULL
    
    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- as.prdlist(lapply(1:length.out, function(i) result), names = conditions)
    #attr(outlist, "pars") <- pars
    
    return(outlist)
    
  }
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- mycondition
  class(outfn) <- c("prdfn", "fn") 
  return(outfn)
  
}

#' Observation function
#'
#' @description An observation function is a function is that is concatenated
#' with a prediction function via \link{prodfn} to yield a new prediction function,
#' see \link{prdfn}. Observation functions are generated, e.g. by \link{Y}. Handling
#' of the conditions is then organized by the \code{obsfn} object.
#' @param X2Y the low-level observation function generated e.g. by \link{Y}.
#' @param parameters character vector with parameter names
#' @param condition character, the condition name
#' @details Observation functions can be "added" by the "+" operator, see \link{sumfn}. Thereby,
#' observations for different conditions are merged or, overwritten. Observation functions can
#' also be concatenated with other functions, e.g. observation functions (\link{obsfn}) or 
#' prediction functions (\link{prdfn}) by the "*" operator, see \link{prodfn}.
#' @return Object of class \code{obsfn}, i.e. a function \code{x(times, pars, fixed, deriv, conditions)}
#' which returns a \link{prdlist}.
#' @export
obsfn <- function(X2Y, parameters = NULL, condition = NULL) {
  
  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- X2Y
  names(mappings) <- condition
  
  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = mycondition) {
   
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
    out <- arglist[[1]]
    pars <- arglist[[2]]
    
     
    overlap <- intersect(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap
    
    if (is.null(overlap) | length(overlap) > 0)
      result <- X2Y(out = out, pars = pars)
    else
      result <- NULL
    
    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- as.prdlist(lapply(1:length.out, function(i) result), names = conditions)
    #attr(outlist, "pars") <- pars
    
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
#' @return Object of class \code{datalist}.
#' @export
datalist <- function(...) {
  mylist <- list(...)
  mynames <- names(mylist)
  if (is.null(mynames)) mynames <- as.character(1:length(mylist))
  as.datalist(mylist, mynames)
}


## Objective classes ---------------------------------------------------------


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

#' Direct sum of objective functions
#' 
#' @param x1 function of class \code{objfn}
#' @param x2 function of class \code{objfn}
#' @details The objective functions are evaluated and their results as added. Sometimes,
#' the evaluation of an objective function depends on results that have been computed
#' internally in a preceding objective function. Therefore, environments are forwarded
#' and all evaluations take place in the same environment. The first objective function
#' in a sum of functions generates a new environment.
#' @return Object of class \code{objfn}.
#' @seealso \link{normL2}, \link{constraintL2}, \link{priorL2}, \link{datapointL2}
#' @aliases sumobjfn
#' @export
"+.objfn" <- function(x1, x2) {
  
  if (is.null(x1)) return(x2)
  
  # objfn + objfn
  if (inherits(x1, "objfn") & inherits(x2, "objfn")) {
    
    outfn <- function(...) {
      args <- list(...)
      v1 <- do.call(x1, args)
      args$env <- attr(v1, "env")
      v2 <- do.call(x2, args)
      out <- v1 + v2
      attr(out, "env") <- args$env
      return(out)
    }
    
    class(outfn) <- c("objfn", "fn")
    return(outfn)
    
  }
  
  
}

#' Direct sum of functions
#'
#' Used to add prediction function, parameter transformation functions or observation functions.
#' 
#' @param x1 function of class \code{obsfn}, \code{prdfn} or \code{parfn}
#' @param x2 function of class \code{obsfn}, \code{prdfn} or \code{parfn}
#' @details Each prediction function is associated to a number of conditions. Adding functions
#' means merging or overwriting the set of conditions.
#' @return Object of the same class as \code{x1} and \code{x2} which returns results for the
#' union of conditions. 
#' @aliases sumfn
#' @seealso \link{P}, \link{Y}, \link{Xs}
#' @export
"+.fn" <- function(x1, x2) {
  
  if (is.null(x1)) return(x2)
   
  mappings.x1 <- attr(x1, "mappings")
  mappings.x2 <- attr(x2, "mappings")
  
  conditions.x1 <- attr(x1, "conditions")
  conditions.x2 <- attr(x2, "conditions")
  overlap <- intersect(conditions.x1, conditions.x2)
  
  
  if (is.null(names(mappings.x1)) || is.null(names(mappings.x2))) stop("General transformations (NULL names) cannot be coerced.")
  
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    mappings.x1 <- mappings.x1[!conditions.x1 %in% overlap]
    conditions.x1 <- conditions.x1[!conditions.x1 %in% overlap]
  }
  
  conditions.x12 <- c(conditions.x1, conditions.x2)
  mappings <- c(mappings.x1, mappings.x2)
  
  # prdfn + prdfn
  if (inherits(x1, "prdfn") & inherits(x2, "prdfn")) {
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings)) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]
      
      
      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)  
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      #outpars <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](times = times, pars = pars, deriv = deriv)
        #outpars[[C]] <- attr(outlist[[C]], "pars")
        #attr(outlist[[C]], "pars") <- NULL
      }
      
      out <- as.prdlist(outlist)
      #attr(out, "pars") <- outpars
      return(out)
      
    }
    
    class(outfn) <- c("prdfn", "fn")
    
  }
  
  # obsfn + obsfn
  if (inherits(x1, "obsfn") & inherits(x2, "obsfn")) {
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings)) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]
      
      
      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)  
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](out = out, pars = pars)
      }
      
      out <- as.prdlist(outlist)
      return(out)
      
    }
    
    class(outfn) <- c("obsfn", "fn")
    
  }
  
  
  # parfn + parfn
  if (inherits(x1, "parfn") & inherits(x2, "parfn")) {
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings)) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]
      
      
      if (is.null(conditions)) {
        available <- names(mappings)
      } else {
        available <- intersect(names(mappings), conditions)  
      }
      outlist <- structure(vector("list", length(conditions)), names = conditions)
      for (C in available) {
        outlist[[C]] <- mappings[[C]](pars = pars, fixed = fixed, deriv = deriv)
      }
      
      return(outlist)
      
    }
    
    class(outfn) <- c("parfn", "fn")
    
  }
  
    
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- union(attr(x1, "parameters"), attr(x2, "parameters"))
  attr(outfn, "conditions") <- conditions.x12
  
  return(outfn)
  
}

out_conditions <- function(c1, c2) {
  
  if (!is.null(c1)) return(c1)
  if (!is.null(c2)) return(c2)
  return(NULL)
  
}


#' Concatenation of functions
#' 
#' Used to concatenate observation functions, prediction functions and parameter transformation functions.
#' 
#' @param p1 function of class \code{obsfn}, \code{prdfn} or \code{parfn}
#' @param p2 function of class \code{obsfn}, \code{prdfn} or \code{parfn}
#' @return Object of the same class as \code{x1} and \code{x2}.
#' @aliases prodfn
#' @export
"*.fn" <- function(p1, p2) {

  # obsfn * obsfn -> obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "obsfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]
      
      
      step1 <- p2(out = out, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = attr(step1[[i]], "parameters"), fixed = fixed, deriv = deriv, conditions = names(step1)[i])))
      
      
      out <- as.prdlist(step2)
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("obsfn", "fn", "composed")
    
    return(outfn)
    
  }  
  
  
  # obsfn * parfn -> obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "parfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]
      
      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = out, pars = step1[[i]], fixed = fixed, deriv = deriv, conditions = names(step1)[i])))
      
      out <- as.prdlist(step2)
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("obsfn", "fn", "composed")
    
    return(outfn)
    
  }  
  
  
  # obsfn * prdfn -> prdfn
  if (inherits(p1, "obsfn") & inherits(p2, "prdfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]
      
      step1 <- p2(times = times, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = attr(step1[[i]], "parameters"), fixed = fixed, deriv = deriv, conditions = names(step1)[i])))
      
      out <- as.prdlist(step2)
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("prdfn", "fn", "composed")
    
    return(outfn)
    
  }  
  
  
  # prdfn * parfn -> prdfn
  if (inherits(p1, "prdfn") & inherits(p2, "parfn")) {
    
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)
    
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]
      
      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(times = times, pars = step1[[i]], deriv = deriv, conditions = names(step1)[i])))
      
      out <- as.prdlist(step2)
      
      return(out)
      
    }
    
    attr(outfn, "mappings") <- attr(p1, "mappings")
    attr(outfn, "conditions") <- conditions.out
    attr(outfn, "parameters") <- attr(p2, "parameters")
    class(outfn) <- c("prdfn", "fn", "composed")
    
    return(outfn)
    
  }
  
  # parfn * parfn -> parfn
  if (inherits(p1, "parfn") & inherits(p2, "parfn")) {
    
    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)
    
    
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL) {
      
      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]
      
      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(pars = step1[[i]], deriv = deriv, conditions = names(step1)[i])))
      return(step2)
      
    }
    
    attr(outfn, "mappings") <- attr(p2, "mappings")
    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("parfn", "fn", "composed")
    
    return(outfn)
    
  }
  
  
}



#' List, get and set controls for different functions
#' 
#' @description Applies to objects of class \code{objfn},
#' \code{parfn}, \code{prdfn} and \code{obsfn}. Allows to manipulate
#' different arguments that have been set when creating the 
#' objects.
#' @details If called without further arguments, \code{controls(x)} lists the
#' available controls within an object. Calling \code{controls()} with \code{name}
#' and \code{condition} returns the control value. The value can be overwritten. If
#' a list or data.frame ist returned, elements of those can be manipulated by the 
#' \code{$}- or \code{[]}-operator.
#' 
#' @param x function
#' @param ... arguments going to the appropriate S3 methods
#' @return Either a print-out or the values of the control. 
#' 
#' @export
controls <- function(x, ...) {
  UseMethod("controls", x)
}



lscontrols_objfn <- function(x) {
  
  names(environment(x)$controls)
  
}

lscontrols_fn <- function(x, condition = NULL) {
  
  conditions <- attr(x, "conditions")
  mappings <- attr(x, "mappings")
  
  
  for (i in 1:length(mappings)) {
    if (is.null(conditions) || is.null(condition) || conditions[i] %in% condition) {
      cat(conditions[i], ":\n", sep = "")
      print(names(environment(mappings[[i]])$controls))  
    }
  }  
  
}

#' @export
#' @rdname controls
#' @param name character, the name of the control
controls.objfn <- function(x, name = NULL, ...) {
  
  if (is.null(name)) lscontrols_objfn(x) else environment(x)$controls[[name]]
}

#' @export
#' @rdname controls
#' @param condition character, the condition name
controls.fn <- function(x, condition = NULL, name = NULL, ...) {
  
  if (is.null(name)) {
    
    lscontrols_fn(x, condition)
    
  } else {
    
    mappings <- attr(x, "mappings")
    if (is.null(condition)) y <- mappings[[1]] else y <- mappings[[condition]]
    environment(y)$controls[[name]]
    
  }
  
}


#' @export
#' @rdname controls
"controls<-" <- function(x, ..., value) {
  UseMethod("controls<-", x)
}


#' @export
#' @param value the new value
#' @rdname controls
"controls<-.objfn" <- function(x, name, ..., value) {
  environment(x)$controls[[name]] <- value
  return(x)
}

#' @export
#' @rdname controls
"controls<-.fn" <- function(x, condition = NULL, name, ..., value) {
  mappings <- attr(x, "mappings")
  if (is.null(condition)) y <- mappings[[1]] else y <- mappings[[condition]]
  environment(y)$controls[[name]] <- value
  return(x)
}