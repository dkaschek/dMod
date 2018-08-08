## ODE model class -------------------------------------------------------------------


#' Generate the model objects for use in Xs (models with sensitivities)
#'
#' @param f Something that can be converted to \link{eqnvec},
#' e.g. a named character vector with the ODE
#' @param deriv logical, generate sensitivities or not
#' @param forcings Character vector with the names of the forcings
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (character or numeric, time point), "value" (character or numeric, value),
#' "method" (character, either
#' "replace" or "add"). See \link[deSolve]{events}. Events need to be defined here if they contain
#' parameters, like the event time or value. If both, time and value are purely numeric, they
#' can be specified in \code{\link{Xs}()}, too.
#' @param outputs Named character vector for additional output variables.
#' @param fixed Character vector with the names of parameters (initial values and dynamic) for which
#' no sensitivities are required (will speed up the integration).
#' @param estimate Character vector specifying parameters (initial values and dynamic) for which
#' sensitivities are returned. If estimate is specified, it overwrites `fixed`.
#' @param modelname Character, the name of the C file being generated.
#' @param solver Solver for which the equations are prepared.
#' @param gridpoints Integer, the minimum number of time points where the ODE is evaluated internally
#' @param verbose Print compiler output to R command line.
#' @param ... Further arguments being passed to funC.
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object)
#' @export
#' @example inst/examples/odemodel.R
#' @import cOde
#' @importFrom digest digest
odemodel <- function(f, deriv = TRUE, forcings=NULL, events = NULL, outputs = NULL, fixed = NULL, estimate = NULL, modelname = "odemodel", solver = c("deSolve", "Sundials"), gridpoints = NULL, verbose = FALSE, ...) {


  if (is.null(gridpoints)) gridpoints <- 2

  f <- as.eqnvec(f)
  # Add hash to prevent overwriting existing models
  modelname <- paste0(modelname, "_", substr(digest(list(f,forcings,events,outputs,fixed,estimate,solver,gridpoints)),1,8))
  modelname_s <- paste0(modelname, "_s")
  solver <- match.arg(solver)

  func <- cOde::funC(f, forcings = forcings, events = events, outputs = outputs, fixed = fixed, modelname = modelname , solver = solver, nGridpoints = gridpoints, ...)
  extended <- NULL
  if (solver == "Sundials") {
    # Sundials does not need "extended" by itself, but dMod relies on it.
    extended <- func
    attr(extended, "deriv") <- TRUE
    attr(extended, "variables") <- c(attr(extended, "variables"), attr(extended, "variablesSens"))
    attr(extended, "events") <- events
  }

  if (deriv && solver == "deSolve") {

    mystates <- attr(func, "variables")
    myparameters <- attr(func, "parameters")

    if (is.null(estimate) & !is.null(fixed)) {
      mystates <- setdiff(mystates, fixed)
      myparameters <- setdiff(myparameters, fixed)
    }

    if (!is.null(estimate)) {
      mystates <- intersect(mystates, estimate)
      myparameters <- intersect(myparameters, estimate)
    }

    s <- sensitivitiesSymb(f,
                           states = mystates,
                           parameters = myparameters,
                           inputs = attr(func, "forcings"),
                           events = attr(func, "events"),
                           reduce = TRUE)
    fs <- c(f, s)
    outputs <- c(attr(s, "outputs"), attr(func, "outputs"))
    events <- rbind(attr(s, "events"), attr(func, "events"))

    extended <- cOde::funC(fs, forcings = forcings, modelname = modelname_s, solver = solver, nGridpoints = gridpoints, events = events, outputs = outputs, ...)
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
#' @seealso \link{eqnlist}
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
#' @example inst/examples/eqnlist.R
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
#' @return object of class \code{parfn}, i.e. a function \code{p(..., fixed, deriv,
#'  conditions, env)}. The argument \code{pars} should be passed via the \code{...}
#'  argument.
#'
#'  Contains attributes "mappings", a list of \code{p2p}
#' functions, "parameters", the union of parameters acceted by the mappings and
#' "conditions", the total set of conditions.
#' @seealso \link{sumfn}, \link{P}
#' @example inst/examples/prediction.R
#' @export
parfn <- function(p2p, parameters = NULL, condition = NULL) {

  force(condition)
  mappings <- list()
  mappings[[1]] <- p2p
  names(mappings) <- condition

  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition, env = NULL) {


    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pars <- arglist[[1]]

    overlap <- test_conditions(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap

    if (is.null(overlap)) conditions <- union(condition, conditions)

    if (is.null(overlap) | length(overlap) > 0)
      result <- p2p(pars = pars, fixed = fixed, deriv = deriv)
    else
      result <- NULL

    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- structure(vector("list", length.out), names = conditions)

    if (is.null(condition)) available <- 1:length.out else available <- match(condition, conditions)
    for (C in available[!is.na(available)]) outlist[[C]] <- result


    return(outlist)

  }
  attr(outfn, "mappings") <- mappings
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- condition
  class(outfn) <- c("parfn", "fn")
  return(outfn)


}

#' Generate a parameter frame
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
#' @details Parameter frames can be subsetted either by \code{[ , ]} or by \code{subset}. If
#' \code{[ , index]} is used, the names of the removed columns will also be removed from
#' the corresponding attributes, i.e. metanames, obj.attributes and parameters.
#' @example inst/examples/parlist.R
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
#' @example inst/examples/parlist.R
#' @seealso \link{load.parlist}, \link{plot.parlist}
parlist <- function(...) {

  mylist <- list(...)
  return(as.parlist(mylist))

}



#' Parameter vector
#'
#' @description A parameter vector is a named numeric vector (the parameter values)
#' together with a "deriv" attribute
#' (the Jacobian of a parameter transformation by which the parameter vector was generated).
#' @param ... objects to be concatenated
#' @param deriv matrix with rownames (according to names of \code{...}) and colnames
#' according to the names of the parameter by which the parameter vector was generated.
#' @return An object of class \code{parvec}, i.e. a named numeric vector with attribute "deriv".
#' @example inst/examples/parvec.R
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
#' @description A prediction function is a function \code{x(..., fixed, deriv, conditions)}.
#' Prediction functions are generated by \link{Xs}, \link{Xf} or \link{Xd}. For an example
#' see the last one.
#'
#' @param P2X transformation function as being produced by \link{Xs}.
#' @param parameters character vector with parameter names
#' @param condition character, the condition name
#' @details Prediction functions can be "added" by the "+" operator, see \link{sumfn}. Thereby,
#' predictions for different conditions are merged or overwritten. Prediction functions can
#' also be concatenated with other functions, e.g. observation functions (\link{obsfn}) or
#' parameter transformation functions (\link{parfn}) by the "*" operator, see \link{prodfn}.
#' @return Object of class \code{prdfn}, i.e. a function \code{x(..., fixed, deriv, conditions, env)}
#' which returns a \link{prdlist}. The arguments \code{times} and
#' \code{pars} (parameter values) should be passed via the \code{...} argument, in this order.
#' @example inst/examples/prediction.R
#' @export
prdfn <- function(P2X, parameters = NULL, condition = NULL) {

  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- P2X
  names(mappings) <- condition

  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = mycondition, env = NULL) {

    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
    times <- arglist[[1]]
    pars <- arglist[[2]]

    # yields derivatives for all parameters in pars but not in fixed
    pars <- c(as.parvec(pars[setdiff(names(pars), names(fixed))]),
              fixed)


    overlap <- test_conditions(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap



    if (is.null(overlap)) conditions <- union(condition, conditions)


    if (is.null(overlap) | length(overlap) > 0)
      result <- P2X(times = times, pars = pars, deriv = deriv)

    else
      result <- NULL

    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- structure(vector("list", length.out), names = conditions)

    if (is.null(condition)) available <- 1:length.out else available <- match(condition, conditions)
    for (C in available[!is.na(available)]) outlist[[C]] <- result
    outlist <- as.prdlist(outlist)

    #length.out <- max(c(1, length(conditions)))
    #outlist <- as.prdlist(lapply(1:length.out, function(i) result), names = conditions)
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
#' see \link{prdfn}. Observation functions are generated by \link{Y}. Handling
#' of the conditions is then organized by the \code{obsfn} object.
#' @param X2Y the low-level observation function generated e.g. by \link{Y}.
#' @param parameters character vector with parameter names
#' @param condition character, the condition name
#' @details Observation functions can be "added" by the "+" operator, see \link{sumfn}. Thereby,
#' observations for different conditions are merged or, overwritten. Observation functions can
#' also be concatenated with other functions, e.g. observation functions (\link{obsfn}) or
#' prediction functions (\link{prdfn}) by the "*" operator, see \link{prodfn}.
#' @return Object of class \code{obsfn}, i.e. a function \code{x(..., fixed, deriv, conditions, env)}
#' which returns a \link{prdlist}. The arguments \code{out} (prediction) and \code{pars} (parameter values)
#' should be passed via the \code{...} argument.
#' @example inst/examples/prediction.R
#' @export
obsfn <- function(X2Y, parameters = NULL, condition = NULL) {

  mycondition <- condition
  mappings <- list()
  mappings[[1]] <- X2Y
  names(mappings) <- condition

  outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = mycondition, env = NULL) {

    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
    out <- arglist[[1]]
    pars <- arglist[[2]]

    # yields derivatives for all parameters in pars but not in fixed
    pars <- c(as.parvec(pars[setdiff(names(pars), names(fixed))]),
              fixed)

    overlap <- test_conditions(conditions, condition)
    # NULL if at least one argument is NULL
    # character(0) if no overlap
    # character if overlap

    if (is.null(overlap)) conditions <- union(condition, conditions)
    if (is.null(overlap) | length(overlap) > 0)
      result <- X2Y(out = out, pars = pars)
    else
      result <- NULL

    # Initialize output object
    length.out <- max(c(1, length(conditions)))
    outlist <- structure(vector("list", length.out), names = conditions)

    if (is.null(condition)) available <- 1:length.out else available <- match(condition, conditions)
    for (C in available[!is.na(available)]) outlist[[C]] <- result
    outlist <- as.prdlist(outlist)

    #length.out <- max(c(1, length(conditions)))
    #outlist <- as.prdlist(lapply(1:length.out, function(i) result), names = conditions)
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
#' parameter vector of inner parameters to produce the prediction frame.
#'
#' Prediction frames are usually the constituents of prediction lists (\link{prdlist}). They are
#' produced by \link{Xs}, \link{Xd} or \link{Xf}. When you define your own prediction functions,
#' see \code{P2X} in \link{prdfn}, the result should be returned as a prediction frame.
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
#' @details Datalists can be plotted, see \link{plotData} and merged, see \link{sumdatalist}.
#' They are the basic structure when combining model prediction and data via the \link{normL2}
#' objective function.
#'
#' Datalists carry the attribute \code{condition.grid} which contains additional information about different
#' conditions, such as dosing information for the experiment. It can be conveniently accessed by the \link{covariates}-function.
#' Reassigning names to a datalist also renames the rows of the \code{condition.grid}.
#' @param ... data.frame objects to be coerced into a list and additional arguments
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
#' @description An objective list contains an objective value, a gradient, and a Hessian matrix.
#'
#' Objective lists can contain additional numeric attributes that are preserved or
#' combined with the corresponding attributes of another objective list when
#' both are added by the "+" operator, see \link{sumobjlist}.
#'
#' Objective lists are returned by objective functions as being generated
#' by \link{normL2}, \link{constraintL2}, \link{priorL2} and \link{datapointL2}.
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
#' @param deriv.err matrix of the derivatives of the error model.
#' @return An object of class \code{objframe}, i.e. a data frame with attribute "deriv".
#' @export
objframe <- function(mydata, deriv = NULL, deriv.err = NULL) {

  # Check column names
  mydata <- as.data.frame(mydata)
  correct.names <- c("time", "name", "value", "prediction",
                     "sigma", "residual", "weighted.residual", "bloq")

  ok <- all(correct.names %in% names(mydata))
  if (!ok) stop("mydata does not have required names")

  out <- mydata[, correct.names]
  attr(out, "deriv") <- deriv
  attr(out, "deriv.err") <- deriv.err
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
#' @example inst/examples/objective.R
#' @export
"+.objfn" <- function(x1, x2) {

  if (is.null(x1)) return(x2)

  conditions.x1 <- attr(x1, "conditions")
  conditions.x2 <- attr(x2, "conditions")
  conditions12 <- union(conditions.x1, conditions.x2)

  parameters.x1 <- attr(x1, "parameters")
  parameters.x2 <- attr(x2, "parameters")
  parameters12 <- union(parameters.x1, parameters.x2)

  modelname.x1 <- attr(x1, "modelname")
  modelname.x2 <- attr(x2, "modelname")
  modelname12 <- union(modelname.x1, modelname.x2)


  # objfn + objfn
  if (inherits(x1, "objfn") & inherits(x2, "objfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = conditions12, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      # 1. If conditions.xi is null, always evaluate xi, but only once
      # 2. If not null, evaluate at intersection with conditions
      # 3. If not null & intersection is empty, don't evaluate xi at all
      v1 <- v2 <- NULL
      if (is.null(conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions.x1, env = env)
      } else if (any(conditions %in% conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, conditions = intersect(conditions, conditions.x1), env = env)
      }

      if (is.null(conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions.x2, env = env)
      } else if (any(conditions %in% conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = intersect(conditions, conditions.x2), env = attr(v1, "env"))
      }

      out <- v1 + v2
      attr(out, "env") <- attr(v1, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    return(outfn)

  }


}


#' Multiplication of objective functions with scalars
#'
#' @description The \code{\%.*\%} operator allows to multiply objects of class objlist or objfn with
#' a scalar.
#'
#' @param x1 object of class objfn or objlist.
#' @param x2 numeric of length one.
#' @return An objective function or objlist object.
#'
#' @export
"%.*%" <- function(x1, x2) {

  if (inherits(x2, "objlist")) {

    out <- lapply(x2, function(x) {
      x1*x
    })
    # Multiply attributes
    out2.attributes <- attributes(x2)[sapply(attributes(x2), is.numeric)]
    attr.names <- names(out2.attributes)
    out.attributes <- lapply(attr.names, function(n) {
      x1*attr(x2, n)
    })
    attributes(out) <- attributes(x2)
    attributes(out)[attr.names] <- out.attributes

    return(out)


  } else if (inherits(x2, "objfn")) {

    conditions12 <- attr(x2, "conditions")
    parameters12 <- attr(x2, "parameters")
    modelname12 <- attr(x2, "modelname")
    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = conditions12, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      v1 <- x1
      v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions, env = attr(v1, "env"))

      out <- v1 %.*% v2
      attr(out, "env") <- attr(v2, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    return(outfn)

  } else {

    x1*x2

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
#' @example inst/examples/prediction.R
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

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings), env = NULL) {

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

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings), env = NULL) {

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

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = names(mappings), env = NULL) {

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

#' Direct sum of datasets
#'
#' Used to merge datasets with overlapping conditions.
#'
#' @param data1 dataset of class \code{datalist}
#' @param data2 dataset of class \code{datalist}
#' @details Each data list contains data frames for a number of conditions.
#' The direct sum of datalist is meant as merging the two data lists and
#' returning the overarching datalist.
#' @return Object of class \code{datalist} for the
#' union of conditions.
#' @aliases sumdatalist
#' @example inst/examples/sumdatalist.R
#' @export
"+.datalist" <- function(data1, data2) {

  overlap <- names(data2)[names(data2) %in% names(data1)]
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    data1 <- data1[!names(data1) %in% names(data2)]
  }

  conditions <- union(names(data1), names(data2))
  data <- lapply(conditions, function(C) rbind(data1[[C]], data2[[C]]))
  names(data) <- conditions

  grid1 <- attr(data1, "condition.grid")
  grid2 <- attr(data2, "condition.grid")

  grid <- combine(grid1, grid2)




  if (is.data.frame(grid)) grid <- grid[!duplicated(rownames(grid)), , drop = FALSE]

  out <- as.datalist(data)
  attr(out, "condition.grid") <- grid

  return(out)
}

out_conditions <- function(c1, c2) {

  if (!is.null(c1)) return(c1)
  if (!is.null(c2)) return(c2)
  return(NULL)

}

test_conditions <- function(c1, c2) {
  if (is.null(c1)) return(NULL)
  if (is.null(c2)) return(NULL)
  return(intersect(c1, c2))
}

#' Concatenation of functions
#'
#' Used to concatenate observation functions, prediction functions and parameter transformation functions.
#'
#' @param p1 function of class \code{obsfn}, \code{prdfn}, \code{parfn} or \code{idfn}
#' @param p2 function of class \code{obsfn}, \code{prdfn}, \code{parfn} or \code{idfn}
#' @return Object of the same class as \code{x1} and \code{x2}.
#' @aliases prodfn
#' @example inst/examples/prediction.R
#' @export
"*.fn" <- function(p1, p2) {

  # obsfn * obsfn -> obsfn
  if (inherits(p1, "obsfn") & inherits(p2, "obsfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]


      step1 <- p2(out = out, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = attr(step1[[i]], "parameters"), fixed = fixed, deriv = deriv, conditions = names(step1)[i])))


      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for observation function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(out, pars) {
        outfn(out = out, pars = pars, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])



      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

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

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("out", "pars"))]
      out <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = out, pars = step1[[i]], fixed = fixed, deriv = deriv, conditions = names(step1)[i])))

      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for observation function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(out, pars) {
        outfn(out = out, pars = pars, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

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

    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(times = times, pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(out = step1[[i]], pars = attr(step1[[i]], "parameters"), fixed = fixed, deriv = deriv, conditions = names(step1)[i])))

      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for prediction function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(times, pars, deriv = TRUE) {
        outfn(times = times, pars = pars, deriv = deriv, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings

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


    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("times", "pars"))]
      times <- arglist[[1]]
      pars <- arglist[[2]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(times = times, pars = step1[[i]], deriv = deriv, conditions = names(step1)[i])))

      out <- as.prdlist(step2)

      return(out)

    }

    # Generate mappings for prediction function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(times, pars, deriv = TRUE) {
        outfn(times = times, pars = pars, deriv = deriv, conditions = conditions.out[i])[[1]]
      }
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings


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


    outfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- do.call(c, lapply(1:length(step1), function(i) p1(pars = step1[[i]], fixed = fixed, deriv = deriv, conditions = names(step1)[i])))
      return(step2)

    }


    # Generate mappings for parameters function
    l <- max(c(1, length(conditions.out)))
    mappings <- lapply(1:l, function(i) {
      mapping <- function(pars, fixed = NULL, deriv = TRUE) {
        outfn(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions.out[i])[[1]]
      }
      m1 <- modelname(p1, conditions = conditions.p1[i])
      m2 <- modelname(p2, conditions = conditions.p2[i])
      attr(mapping, "modelname") <- union(m1, m2)
      attr(mapping, "parameters") <- getParameters(p2, conditions = conditions.out[i])

      return(mapping)
    })
    names(mappings) <- conditions.out
    attr(outfn, "mappings") <- mappings


    attr(outfn, "parameters") <- attr(p2, "parameters")
    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("parfn", "fn", "composed")

    return(outfn)

  }


  # objfn * parfn -> objfn
  if (inherits(p1, "objfn") & inherits(p2, "parfn")) {

    conditions.p1 <- attr(p1, "conditions")
    conditions.p2 <- attr(p2, "conditions")
    conditions.out <- out_conditions(conditions.p1, conditions.p2)

    outfn <- function(...,  fixed = NULL, deriv=TRUE, conditions = NULL, env = NULL) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, "pars")]
      pars <- arglist[[1]]

      step1 <- p2(pars = pars, fixed = fixed, deriv = deriv, conditions = conditions)
      step2 <- Reduce("+", lapply(1:length(step1), function(i) p1(pars = step1[[i]], fixed = NULL, deriv = deriv, conditions = names(step1)[i], env = env)))
      return(step2)


    }

    attr(outfn, "conditions") <- conditions.out
    class(outfn) <- c("objfn", "fn", "composed")

    return(outfn)


  }

  # idfn * fn -> fn
  if (inherits(p1, "idfn")) {
    return(p2)
  }

  # fn * idfn -> fn
  if (inherits(p2, "idfn")) {
    return(p1)
  }

}

## General purpose functions for different dMod classes ------------------------------

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
#' @examples
#' ## parfn with condition
#' p <- P(eqnvec(x = "-a*x"), method = "implicit", condition = "C1")
#' controls(p)
#' controls(p, "C1", "keep.root")
#' controls(p, "C1", "keep.root") <- FALSE
#'
#' ## obsfn with NULL condition
#' g <- Y(g = eqnvec(y = "s*x"), f = NULL, states = "x", parameters = "s")
#' controls(g)
#' controls(g, NULL, "attach.input")
#' controls(g, NULL, "attach.input") <- FALSE
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


#' Extract the derivatives of an object
#'
#' @param x object from which the derivatives should be extracted
#' @param ... additional arguments (not used right now)
#' @return The derivatives in a format that depends on the class of \code{x}.
#' This is
#' \code{parvec -> matrix},
#' \code{prdframe -> prdframe},
#' \code{prdlist -> prdlist},
#' \code{objlist -> named numeric}.
#' @export
getDerivs <- function(x, ...) {
  UseMethod("getDerivs", x)
}

#' @export
#' @rdname getDerivs
getDerivs.parvec <- function(x, ...) {

  attr(x, "deriv")

}

#' @export
#' @rdname getDerivs
getDerivs.prdframe <- function(x, ...) {

  prdframe(prediction = attr(x, "deriv"), parameters = attr(x, "parameters"))

}

#' @export
#' @rdname getDerivs
getDerivs.prdlist <- function(x, ...) {

  as.prdlist(
    lapply(x, function(myx) {
      getDerivs(myx)
    }),
    names = names(x)
  )

}

#' @export
#' @rdname getDerivs
getDerivs.list <- function(x, ...) {

  lapply(x, function(myx) getDerivs(myx))

}


#' @export
#' @rdname getDerivs
getDerivs.objlist <- function(x, ...) {

  x$gradient

}



getEquations <- function(x, ...) {

  UseMethod("getEquations", x)

}



#' Extract the parameters of an object
#'
#' @param ... objects from which the parameters should be extracted
#' @param conditions character vector specifying the conditions to
#' which \code{getParameters} is restricted
#' @return The parameters in a format that depends on the class of \code{x}.
#' @export
getParameters <- function(..., conditions = NULL) {


  Reduce("union", lapply(list(...), function(x) {
    UseMethod("getParameters", x)
  }))


}



#' @export
#' @rdname getParameters
#' @param x object from which the parameters are extracted
getParameters.odemodel <- function(x, conditions = NULL) {

  parameters <- c(
    attr(x$func, "variables"),
    attr(x$func, "parameters")
  )

  return(parameters)

}


#' @export
#' @rdname getParameters
getParameters.fn <- function(x, conditions = NULL) {

  if (is.null(conditions)) {
    parameters <- attr(x, "parameters")
  } else {
    mappings <- attr(x, "mappings")
    mappings <- mappings[intersect(names(mappings), conditions)]
    parameters <- Reduce("union",
                         lapply(mappings, function(m) attr(m, "parameters"))
    )
  }

  return(parameters)

}
#' @export
#' @rdname getParameters
getParameters.parvec <- function(x, conditions = NULL) {

  names(x)

}

#' @export
#' @rdname getParameters
getParameters.prdframe <- function(x, conditions = NULL) {

  attr(x, "parameters")

}

#' @export
#' @rdname getParameters
getParameters.prdlist <- function(x, conditions = NULL) {

  select <- 1:length(x)
  if (!is.null(conditions)) select <- intersect(names(x), conditions)
  lapply(x[select], function(myx) getParameters(myx))

}

#' @export
#' @rdname getParameters
getParameters.eqnlist <- function(x) {
  unique(c(getSymbols(x$states), getSymbols(x$rates), getSymbols(x$volumes)))
}


#' Extract the conditions of an object
#'
#' @param x object from which the conditions should be extracted
#' @param ... additional arguments (not used right now)
#' @return The conditions in a format that depends on the class of \code{x}.
#' @export
getConditions <- function(x, ...) {
  UseMethod("getConditions", x)
}


#' @export
#' @rdname getConditions
getConditions.list <- function(x, ...) {

  names(x)

}


#' @export
#' @rdname getConditions
getConditions.fn <- function(x, ...) {

  attr(x, "conditions")

}

#' Get and set modelname
#'
#' @description The modelname attribute refers to the name of a C file associated with
#' a dMod function object like prediction-, parameter transformation- or
#' objective functions.
#'
#' @param ... objects of type \code{prdfn}, \code{parfn}, \code{objfn}
#' @param conditions character vector of conditions
#' @return character vector of model names, corresponding to C files
#' in the local directory.
#'
#' @export
modelname <- function(..., conditions = NULL) {

  Reduce("union", lapply(list(...), mname, conditions = conditions))

}

#' Get modelname from single object (used internally)
#'
#' @param x dMod object
#' @param conditions character vector of conditions
#' @export
mname <- function(x, conditions = NULL) {
  UseMethod("mname", x)
}

#' @export
#' @rdname mname
mname.NULL <- function(x, conditions = NULL) NULL

#' @export
#' @rdname mname
mname.character <- function(x, conditions = NULL) {

  mname(get(x), conditions = conditions)

}

#' @export
#' @rdname mname
mname.objfn <- function(x, conditions = NULL) {

  attr(x, "modelname")

}

#' @export
#' @rdname mname
mname.fn <- function(x, conditions = NULL) {

  mappings <- attr(x, "mappings")
  select <- 1:length(mappings)
  if (!is.null(conditions)) select <- intersect(names(mappings), conditions)
  modelnames <- Reduce("union",
                       lapply(mappings[select], function(m) attr(m, "modelname"))
  )

  return(modelnames)

}

#' @export
#' @rdname modelname
#' @param x dMod object for which the model name should be set
#' @param value character, the new modelname (does not change the C file)
"modelname<-" <- function(x, ..., value) {
  UseMethod("modelname<-", x)
}

#' @export
#' @rdname modelname
"modelname<-.fn" <- function(x, conditions = NULL, ..., value) {

  mappings <- attr(x, "mappings")
  select <- 1:length(mappings)
  if (!is.null(conditions)) select <- intersect(names(mappings), conditions)
  #if (length(value) > 1 && length(value) != length(mappings[select]))
  #  stop("Length of modelname vector should be either 1 or equal to the number of conditions.")
  if (length(value) == 1) {
    value <- rep(value, length.out = length(mappings[select]))
    if (!is.null(conditions)) names(value) <- conditions
  }


  for (i in select) {
    attr(attr(x, "mappings")[[i]], "modelname") <- value[i]
    if (inherits(x, "prdfn")) {
      attr(environment(attr(x, "mappings")[[i]])[["extended"]], "modelname") <- value[i]
      attr(environment(attr(x, "mappings")[[i]])[["func"]], "modelname") <- value[i]
    }
  }


  return(x)

}

#' @export
#' @rdname modelname
"modelname<-.objfn" <- function(x, conditions = NULL, ..., value) {
  attr(x, "modelname") <- value
  return(x)
}




#' Extract the equations of an object
#'
#' @param x object from which the equations should be extracted
#' @param conditions character or numeric vector specifying the conditions to
#' which \code{getEquations} is restricted. If \code{conditions} has length one,
#' the result is not returned as a list.
#' @return The equations as list of \code{eqnvec} objects.
#' @export
getEquations <- function(x, conditions = NULL) {

    UseMethod("getEquations", x)

}



#' @export
#' @rdname getEquations
getEquations.odemodel <- function(x, conditions = NULL) {

  attr(x$func, "equations")

}



#' @export
#' @rdname getEquations
getEquations.prdfn <- function(x, conditions = NULL) {

  mappings <- attr(x, "mappings")

  if (is.null(conditions)) {
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    return(equations)
  }

  if (!is.null(conditions)) {
    mappings <- mappings[conditions]
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    if (length(equations) == 1) {
      return(equations[[1]])
    } else {
      return(equations)
    }
  }

}


#' @export
#' @rdname getEquations
getEquations.fn <- function(x, conditions = NULL) {

  mappings <- attr(x, "mappings")

  if (is.null(conditions)) {
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    return(equations)
  }

  if (!is.null(conditions)) {
    mappings <- mappings[conditions]
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    if (length(equations) == 1) {
      return(equations[[1]])
    } else {
      return(equations)
    }
  }

}

#' Extract the observables of an object
#'
#' @param x object from which the equations should be extracted
#' @param ... not used
#' @return The equations as a character.
#' @export
getObservables <- function(x, ...) {
  UseMethod("getObservables", x)
}

