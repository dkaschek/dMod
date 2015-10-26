## Equation classes -------------------------------------------------------

#' Generate equation vector object
#' 
#' @description The eqnvec object stores explicit algebraic equations, like the
#' right-hand sides of an ODE, observation functions or parameter transformations
#' as named character vectors.
#' @param equations (named) character of symbolic mathematical expressions, 
#' the right-hand sides of the equations
#' @param names character, the left-hand sides of the equation
#' @return object of class \code{eqnvec}, basically a named character.
#' @export
eqnvec <- function(equations = NULL, names = NULL) {
  
  if (is.null(equations)) return(c())
  
  if (is.null(names)) names <- names(equations)
  if (is.null(names)) stop("equations need names")
  if (length(names) != length(equations)) stop("Length of names and equations do not coincide")
  try.parse <- try(parse(text = equations), silent = TRUE)
  if (inherits(try.parse, "try-error")) stop("equations cannot be parsed")
  
  out <- structure(equations, names = names)
  class(out) <- "eqnvec"
  
  return(out)
  
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
  class(out) <- "eqnlist"
  
  return(out)
}



## Parameter classes --------------------------------------------------------




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

#' @export
parlist <- function(mylist, myframe) {
  
  # Wie erzeugt man myframe aus mylist?
  
  out <- mylist
  attr(out, "parframe") <- myframe
  
  return(out)
  
}

#' @export
parvec <- function(p, mynames = names(p), deriv = NULL) {
  
  out <- as.numeric(p)
  names(out) <- mynames
  attr(out, "deriv") <- deriv
  class(out) <- c("parvec", "numeric")
  
  return(out)
  
}


## Prediction classes ----------------------------------------------------

#' @export
prdfn <- function(..., pouter = NULL, conditions = "1") {
  
  force(pouter)
  force(conditions)
  
  myexpr <- as.expression(substitute(...))
  
  # Dummy constructor
  is.nullexpr <- deparse(myexpr)[1] == "expression(NULL)"
  if (is.nullexpr) myexpr <- expression({
    out <- matrix(times, ncol = 1)
    colnames(out) <- "time"
    myderivs <- out
    prdframe(out, myderivs, names(pars))
  })
  
  # Prediction function
  myfn <- function(times, pars = pouter, fixed = NULL, deriv = TRUE, ...){
    
    prdlist(
      lapply(conditions, function(condition) {
        eval(myexpr)    
      }),
      conditions
    )
  }
  class(myfn) <- "prdfn"
  attr(myfn, "pouter") <- pouter
  attr(myfn, "parameters") <- names(pouter)
  return(myfn)
  
}


#' @export
prdframe <- function(prediction = NULL, deriv = NULL, sensitivities = NULL, parameters = NULL) {
  
  out <- if (!is.null(prediction)) as.matrix(prediction) else matrix()
  
  attr(out, "deriv") <- deriv
  attr(out, "sensitivities") <- sensitivities
  attr(out, "parameters") <- parameters
  class(out) <- c("prdframe", "matrix")
  
  return(out)
  
}

#' @export
prdlist <- function(mylist = NULL, mynames = names(mylist)) {
  
  if (is.null(mylist)) mylist <- list()
  
  if (length(mynames) != length(mylist)) stop("names argument has wrong length")
  
  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- c("prdlist", "list")
  
  return(mylist)
  
}



## Data classes ----------------------------------------------------------------

#' Generate a datalist object
#' 
#' @description The datalist object stores time-course data in a list of data.frames.
#' The names of the list serve as identifiers, e.g. of an experimental condition, etc.
#' @param mylist list of data.frame, each data.frame is expected to have columns "name" (factor or character),
#' "time" (numeric), "value" (numeric) and "sigma" (numeric).
#' @param mynames character vector of the length of mylist.
#' @return Object of class \code{datalist}.
#' @export
datalist <- function(mylist, mynames = names(mylist)) {
  
  ## Check properties
  is.data.frame <- sapply(mylist, class) == "data.frame"
  if (!all(is.data.frame)) stop("list of data.frame expected")
  
  correct.names <- c("name", "time", "value", "sigma")
  have.correct.names <- sapply(mylist, function(d) all(colnames(d) %in% correct.names))
  if (!all(have.correct.names)) stop(paste("data.frames should have names:", correct.names, collapse = " "))
  
  if (length(mynames) != length(mylist)) stop("names argument has wrong length")
  
  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- "datalist"
  
  return(mylist)
  
}


## Objective classes ---------------------------------------------------------


#' @export
objfn <- function(..., data, x, pouter = NULL, conditions = names(data)) {
  
  force(data)
  force(x)
  force(pouter)
  force(conditions)
  
  myexpr <- as.expression(substitute(...))
  timesD <- sort(c(0, unique(do.call(c, lapply(data, function(d) d$time)))))
  
  myfn <- function(pouter = pouter, fixed = NULL, deriv=TRUE){
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv)
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) wrss(res(data[[cn]], prediction[[cn]])))
    out.data <- Reduce("+", out.data)
    
    # Evaluate user-defined contributions, e.g. priors
    out.user <- eval(myexpr)
    #attributes.user <- attributes(out.user)
    #attributes.user <- attributes.user[!names(attributes.user) %in% c("names", "class")]
    
    # Combine contributions and attach attributes
    out <- out.data + out.user
    attr(out, "data") <- out.data$value
    attr(out, "user") <- out.user$value
    return(out)
    
      
  }
  class(myfn) <- "objfn"
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
  class(out) <- "objlist"
  return(out)
  
}


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
  class(out) <- "objframe"
  
  return(out)
  
  
}
