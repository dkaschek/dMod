
#' Model Predictions
#' 
#' Make a model prediction for times and a parameter frame. The
#' function is a generalization of the standard prediction by a
#' prediction function object in that it allows to pass a parameter
#' frame instead of a single parameter vector.
#' 
#' @param object prediction function
#' @param ... Further arguments goint to the prediction function
#' @param times numeric vector of time points
#' @param pars parameter frame, e.g. output from \link{mstrust} or 
#' \link{profile}
#' @param data data list object. If data is passed, its condition.grid
#' attribute is used to augment the output dataframe by additional 
#' columns. \code{"data"} itself is returned as an attribute.
#' @return A data frame
#' @export
predict.prdfn <- function(object, ..., times, pars, data = NULL) {
  
  
  x <- object
  arglist <- list(...)
  if (any(names(arglist) == "conditions")) {
    C <- arglist[["conditions"]]
    if (!is.null(data)) {
      data <- data[C]
    }
  }
  if (is.null(data)) data <- data.frame()
  condition.grid.data <- attr(data, "condition.grid")
  
  prediction <- do.call(combine, lapply(1:nrow(pars), function(i) {
    
    mypar <- as.parvec(pars, i)
    prediction <- x(times, mypar, deriv = FALSE, ...)
    
    if (is.null(names(prediction))) {
      conditions <- 1
    } else {
      conditions <- names(prediction)
    }
    
    condition.grid <- data.frame(row.names = conditions)
    
    # Augment by parframe metanames and obj.attributes
    mygrid <- pars[i, !colnames(pars) %in% attr(pars, "parameters")]
    mynames <- colnames(mygrid)
    if (length(mynames) > 0) {
      mynames <- paste0(".", mynames)
      colnames(mygrid) <- mynames
      condition.grid <- cbind(condition.grid, mygrid)
    }
    
    # Augment by condition.grid of data
    if (!is.null(condition.grid.data) && ncol(condition.grid.data) > 1) 
      condition.grid <- cbind(condition.grid.data[conditions,], condition.grid)
    
    # Write condition.grid into data
    attr(data, "condition.grid") <- condition.grid

    # Return
    as.data.frame(prediction, data = data)
    
  }))
  
  n <- nrow(prediction)
  
  if (length(data) > 0) {
    attr(data, "condition.grid") <- condition.grid.data
    data <- as.data.frame(data)
    tmp <- combine(prediction, data)
    data <- tmp[-(1:n),]
  }
  
  attr(prediction, "data") <- data
  return(prediction)  
  
  
  
}

#' Generate sample for multi-start fit
#' 
#' @param center named numeric, the center around we sample
#' @param samplefun character, indicating the random number generator,
#' defaults to \code{"rnorm"}.
#' @param fits length of the sample
#' @param ... arguments going to \code{samplefun}
#' @return matrix with the parameter samples
#' @export
mssample <- function(center, samplefun = "rnorm", fits = 20, ...) {
  
  sample.matrix <- do.call(rbind, lapply(1:fits, function(i) 
    do.call(samplefun, c(list(n = length(center), mean = center), list(...)))))
  
  colnames(sample.matrix) <- names(center)
  
  myframe <- parframe(
    cbind(
      index = 1:nrow(sample.matrix),
      as.data.frame(sample.matrix)
    ),
    metanames = "index",
    parameters= colnames(sample.matrix)
  )
  
  return(myframe)

}

#' Reparameterization
#' 
#' @param expr character of the form \code{"lhs ~ rhs"} where \code{rhs}
#' reparameterizes \code{lhs}. Both \code{lhs} and \code{rhs}
#' can contain a number of symbols whose values need to be passed by the \code{...} argument.
#' @param trafo character or equation vector or list thereof. The object where the replacement takes place in
#' @param ... pass symbols as named arguments
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
repar <- function(expr, trafo = NULL, ...) {
 
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
    trafo <- structure(lhs, names = lhs)
  } else if (is.list(trafo)) {
    trafo <- lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  } else if (is.character(trafo)) {
    trafo <- replaceSymbols(lhs, rhs, trafo)
  }
  
  return(trafo)
  
  
}

