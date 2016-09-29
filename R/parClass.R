## Methods for the class parlist -----------------------------------------------

#' Parameter list
#' 
#' @param x list of lists, as returned by \code{trust}
#' @rdname parlist
#' @export
as.parlist <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  } else {
    class(x) <- c("parlist", "list")
    return(x)
  }
}

#' @export
#' @param object a parlist
#' @rdname parlist
summary.parlist <- function(object, ...) {
  
  x <- object
  
  # Statistics
  m_stat <- stat.parlist(x)
  m_error <- sum(m_stat == "error")
  m_converged <- sum(m_stat == "converged")
  m_notConverged <- sum(m_stat == "notconverged")
  m_sumStatus <- sum(m_error, m_converged, m_notConverged)
  m_total <- length(m_stat)
  
  # Best and worst fit
  m_parframe <- as.parframe(x)
  m_order <- order(m_parframe$value)
  m_bestWorst <- m_parframe[c(m_order[1], tail(m_order, 1)),]
  rownames(m_bestWorst) <- c("best", "worst")
  cat("Results of the best and worst fit\n")
  print(m_bestWorst)
  
  cat("\nStatistics of fit outcome",
      "\nFits aborted:       ", m_error,
      "\nFits not converged: ", m_notConverged,
      "\nFits converged:     ", m_converged,
      "\nFits total:         ", m_sumStatus, " [", m_total, "]", sep = "")
}



#' Gather statistics of a fitlist
#' @param x The fitlist
stat.parlist <- function(x) {
  status <- do.call(rbind, lapply(x, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error" || any(is.null(fit)))) {
      return("error")
    } else {
      if (fit$converged) {
        return("converged")
      } else {
        return("notconverged")
      }
    }
  }))
  
  rownames(status) <- 1:length(status)
  colnames(status) <- "fit status"
  
  return(status)
}


#' Coerce object to a parameter frame
#' 
#' @param x object to be coerced
#' @param ... other arguments
#' @return object of class \link{parframe}.
#' @example inst/examples/parlist.R
#' @export
as.parframe <- function(x, ...) {
  UseMethod("as.parframe", x)
}

#' @export
#' @rdname as.parframe
#' @param sort.by character indicating by which colum the returned parameter frame
#' should be sorted. Defaults to \code{"value"}.
as.parframe.parlist <- function(x, sort.by = "value", ...) {
  m_stat <- stat.parlist(x)
  m_metanames <- c("index", "value", "converged", "iterations")
  m_idx <- which("error" != m_stat)
  m_parframe <- do.call(rbind,
                        mapply(function(fit, idx) {
                          data.frame(
                            index = idx,
                            value = fit$value,
                            converged = fit$converged,
                            iterations = fit$iterations,
                            as.data.frame(as.list(fit$argument))
                          )
                        }, fit = x[m_idx], idx = m_idx, SIMPLIFY = FALSE))
  # Sort by value
  m_parframe <- m_parframe[order(m_parframe[sort.by]),]
  
  parframe(m_parframe, parameters = names(x[[m_idx[1]]]$argument), metanames = m_metanames)
  
  
}


## Methods for the class parframe -----------------------------------------------

#' @export
#' @param i row index in any format
#' @param j column index in any format
#' @param drop logical. If TRUE the result is coerced to the lowest possible dimension
#' @rdname parframe
"[.parframe" <- function(x, i = NULL, j = NULL, drop = FALSE){
  
  metanames <- attr(x, "metanames")
  obj.attributes <- attr(x, "obj.attributes")
  parameters <- attr(x, "parameters")
  
  out <- as.data.frame(x)
  #out <- as.data.frame(unclass(x))
  if (!is.null(i)) out <- out[i, ]
  if (!is.null(j)) out <- out[, j, drop = drop]
  
  if (drop) return(out)
  
  metanames <- intersect(metanames, colnames(out))
  obj.attributes <- intersect(obj.attributes, colnames(out))
  parameters <- intersect(parameters, colnames(out))
  
  parframe(out, parameters = parameters, metanames = metanames, obj.attributes = obj.attributes)
  
}


#' @export
#' @param ... additional arguments
#' @rdname parframe
subset.parframe <- function(x, ...) {
  
  x[with(as.list(x), ...), ]
  
}




## Methods for the class parvec ------------------------------------------------

#' Dispatch as.parvec.
#'
#' @export
#' @rdname parvec
as.parvec <- function(x, ...) {
  UseMethod("as.parvec", x)
}


#' Parameter vector
#' @param x numeric or named numeric, the parameter values
#' @param names optional character vector, the parameter names. Otherwise, names
#' are taken from \code{x}.
#' @rdname parvec
#' @export
as.parvec.numeric <- function(x, names = NULL, deriv = NULL, ...) {
  
  p <- x
  
  out <- as.numeric(p)
  if (is.null(names)) names(out) <- names(p) else names(out) <- names
  if (is.null(deriv)) deriv <- attr(x, "deriv")
  if (is.null(deriv)) {
    deriv <- diag(length(out))
    colnames(deriv) <- rownames(deriv) <- names(out)
  }
  attr(out, "deriv") <- deriv
  class(out) <- c("parvec", "numeric")
  
  return(out)
  
}


#' Pretty printing for a parameter vector
#' 
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @param x object of class \code{parvec}
#' @param ... not used yet.
#' @export
print.parvec <- function(x, ...) {
  
  par <- x
  
  m_parWidth <- max(nchar(names(par)))
  m_names <- names(par)
  m_order <- order(m_names)
  
  msg <- mapply(function(p, n) {
    if (!as.numeric(p) < 0 ) {
      p <- paste0(" ", p)
    }
    paste0(strpad(n, m_parWidth, where = "left"), ": ", p)
  }, p = par[m_order], n = m_names[m_order])
  
  cat(msg, sep = "\n")
}




#' @export
#' @param drop logical, drop empty columns in Jacobian after subsetting. 
#' ATTENTION: Be careful with this option. The default behavior is to keep
#' the columns in the Jacobian. This can lead to unintended results when
#' subsetting the parvec and using it e.g. in another parameter
#' transformation.
#' @rdname parvec
"[.parvec" <- function(x, ..., drop = FALSE) {
  out <- unclass(x)[...]
  deriv <- submatrix(attr(x, "deriv"), rows = ...)
  if (drop) {
    empty.cols <- apply(deriv, 2, function(v) all(v == 0))
    deriv <- submatrix(deriv, cols = !empty.cols)
  }
  as.parvec(out, deriv = deriv)
}

#' @export
#' @rdname parvec
c.parvec <- function(...) {
  
  mylist <- list(...) #lapply(list(...), as.parvec)
  
  n <- unlist(lapply(mylist, function(l) names(l)))
  v <- unlist(lapply(mylist, function(l) as.numeric(l)))
  d <- lapply(mylist, function(l) attr(l, "deriv"))
  
  
  
  if (any(duplicated(n))) stop("Found duplicated names. Parameter vectors cannot be coerced.")
  
  deriv <- Reduce(combine, d)
  n.missing <- setdiff(n, rownames(deriv))
  n.available <- intersect(n, rownames(deriv))
  deriv.missing <- matrix(0, nrow = length(n.missing), ncol = ncol(deriv), 
                          dimnames = list(n.missing, colnames(deriv)))
  
  ## Attention: The expected way of function is that
  ## no columns are attachd for parameters for which no derivatives
  ## were available. This is important for prdfn() and obsfn() to 
  ## work properly with the "fixed" argument.
  deriv <- submatrix(rbind(deriv, deriv.missing), rows = n)
  
  as.parvec(v, names = n, deriv = deriv)
  
}



## Methods for the class parfn--------------------------------------------------

#' Pretty printing parameter transformations
#' 
#' @param x prediction function
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
print.parfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Parameter transformation:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  
}

#' @export
summary.parfn <- function(object, ...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}


