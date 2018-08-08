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
    if (inherits(fit, "try-error") || any(names(fit) == "error") || any(is.null(fit))) {
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


#' Plot a parameter list.
#' 
#' @param x fitlist obtained from mstrust
#' @param ... additional arguments
#' @param path print path of parameters from initials to convergence. For this
#'   option to be TRUE \code{\link{mstrust}} must have had the option
#'   \option{blather}.
#' 
#' @details If path=TRUE:        
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
plot.parlist <- function(x, path = FALSE, ...) {
  
  pl <- x
  
  index <- do.call(rbind, lapply(pl, function(l) l$converged))
  fl <- pl[index]
  if (!path) {
    initPar <- do.call(rbind, lapply(fl, function(l) l$parinit))
    convPar <- do.call(rbind, lapply(fl, function(l) l$argument))
    
    ddata <- data.frame(cbind(matrix(initPar, ncol = 1), matrix(convPar, ncol = 1) ))
    ddata <- cbind(rep(colnames(initPar), each = nrow(initPar)), ddata, 1)
    names(ddata) <- c("parameter","x","y","run")
    
    #plot initial vs converged parameter values
    ggplot(data=ddata)+facet_wrap(~ parameter)+geom_point(aes(x=x,y=y))
  } else {
    if (!any (names(fl[[1]]) == "argpath")){
      stop("No path information in the output of mstrust. Restart mstrust with option blather.")
    }
    parNames <- names(fl[[1]]$parinit)
    
    pathPar <- do.call(rbind, mapply(function(l, idx) {
      mParPath <- as.data.frame(matrix(l$argpath, ncol = 1))
      mParPath <- cbind(rep(parNames,each = nrow(l$argpath), times = 1), rep(1:nrow(l$argpath), length(parNames)), mParPath, as.character(idx))
    }, l = fl, idx = 1:length(fl), SIMPLIFY = FALSE))
    names(pathPar) <- c("parameter", "iteration", "path", "idx")
    ggplot(data=pathPar)+geom_line(aes(x=iteration,y=path,colour=idx))+facet_wrap(~ parameter)
  }
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



#' Concatenate parameter lists
#'
#' @description Fitlists carry an fit index which must be held unique on merging
#' multiple fitlists.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @rdname parlist
#' @export
#' @export c.parlist
c.parlist <- function(...) {
    m_fits <- lapply(list(...), unclass)
    m_fits <- do.call(c, m_fits)
    m_parlist <- mapply(function(fit, idx) {
      if (is.list(fit)) fit$index <- idx
      return(fit)
      }, fit = m_fits, idx = seq_along(m_fits), SIMPLIFY = FALSE)
    
    return(as.parlist(m_parlist))
  }





## Methods for the class parframe ----


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


#' Select a parameter vector from a parameter frame.
#' 
#' @description Obtain a parameter vector from a parameter frame.
#' 
#' @param x A parameter frame, e.g., the output of
#'   \code{\link{as.parframe}}.
#' @param index Integer, the parameter vector with the \code{index}-th lowest
#'   objective value.
#' @param ... not used right now
#'   
#' @details With this command, additional information included in the parameter
#'   frame as the objective value and the convergence state are removed and a
#'   parameter vector is returned. This parameter vector can be used to e.g.,
#'   evaluate an objective function.
#'   
#'   On selection, the parameters in the parameter frame are ordered such, that
#'   the parameter vector with the lowest objective value is at \option{index}
#'   1. Thus, the parameter vector with the \option{index}-th lowest objective
#'   value is easily obtained.
#'   
#' @return The parameter vector with the \option{index}-th lowest objective
#'   value.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
as.parvec.parframe <- function(x, index = 1, ...) {
  parframe <- x
  m_order <- 1:nrow(x)
  metanames <- attr(parframe, "metanames")
  if ("value" %in% metanames) m_order <- order(parframe$value)
  best <- as.parvec(unlist(parframe[m_order[index], attr(parframe, "parameters")]))
  if ("converged" %in% metanames && !parframe[m_order[index],]$converged) {
    warning("Parameter vector of an unconverged fit is selected.", call. = FALSE)
    }
  return(best)
}




#' @export
#' @rdname plotPars
plotPars.parframe <- function(x, tol = 1, ...){
  
  if (!missing(...)) x <- subset(x, ...)
  
  jumps <- stepDetect(x$value, tol)
  jump.index <- approx(jumps, jumps, xout = 1:length(x$value), method = "constant", rule = 2)$y
  
  #values <- round(x$value/tol)
  #unique.values <- unique(values)
  #jumps <- which(!duplicated(values))
  #jump.index <- jumps[match(values, unique.values)]
  x$index <- as.factor(jump.index)
  
  myparframe <- x
  parNames <- attr(myparframe,"parameters")
  parOut <- wide2long.data.frame(out = ((myparframe[, c("index", "value", parNames)])) , keep = 1:2)
  names(parOut) <- c("index", "value", "name", "parvalue")
  plot <- ggplot2::ggplot(parOut, aes(x = name, y = parvalue, color = index)) + geom_boxplot(outlier.alpha = 0) + theme_dMod() + scale_color_dMod() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  attr(plot, "data") <- parOut
  
  return(plot)
  
}


#' @export
#' @rdname plotValues
plotValues.parframe <- function(x, tol = 1, ...) {
  
  if (!missing(...)) x <- subset(x, ...)
  
  jumps <- stepDetect(x$value, tol)
  y.jumps <- seq(max(x$value), min(x$value), length.out = length(jumps))
  
  
  pars <- x
  pars <- pars[order(pars$value),]
  pars[["index"]] <-  1:nrow(pars)
  
  P <- ggplot2::ggplot(pars, aes(x = index, y = value, pch = converged, color = iterations)) + 
    geom_vline(xintercept = jumps, lty = 2) +
    geom_point() + 
    annotate("text", x = jumps + 1, y = y.jumps, label = jumps, hjust = 0, color = "red", size = 3) +
    xlab("index") + ylab("value") + theme_dMod()
  
  attr(P, "data") <- pars
  attr(P, "jumps") <- jumps
  
  return(P)
  
}



#' @export
#' @rdname plotProfile
plotProfile.parframe <- function(profs, ..., maxvalue = 5, parlist = NULL) {
  
  if("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])
    
    
    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }
    
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue
      
      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)
      
      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }
      
      return(sub)
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))
  
  data.zero <- subset(data, is.zero)
  
  threshold <- c(1, 2.7, 3.84)
  
  data <- droplevels.data.frame(subset(data, ...))
  
  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") + 
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")
  
  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){ 
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]  
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)
    
    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)
  
}


#' @export
#' @rdname plotProfile
plotProfile.list <- function(profs, ..., maxvalue = 5, parlist = NULL) {
  
  if("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])
    
    
    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }
    
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue
      
      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)
      
      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }
      
      return(sub)
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))
  
  data.zero <- subset(data, is.zero)
  
  threshold <- c(1, 2.7, 3.84)
  
  data <- droplevels.data.frame(subset(data, ...))
  
  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") + 
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")
  
  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){ 
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]  
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)
    
    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)
  
}



#' @export
#' @rdname parframe
is.parframe <- function(x) {
  "parframe" %in% class(x)
}

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

#' Extract those lines of a parameter frame with unique elements in the value column
#' @param x parameter frame
#' @param incomparables not used. Argument exists for compatibility with S3 generic.
#' @param tol tolerance to decide when values are assumed to be equal, see \code{\link{plotValues}()}.
#' @param ... additional arguments being passed to \code{\link{plotValues}()}, e.g. for subsetting.
#' @return A subset of the parameter frame \code{x}.
#' @export
unique.parframe <- function(x, incomparables = FALSE, tol = 1, ...) {
  
  
  jumps <- attr(plotValues(x = x, tol = tol, ...), "jumps")
  x[jumps, ]
  
  
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
  
  # myclass <- class(...)
  # if (inherits(myclass, "character")) {
  #   select.name <- Reduce("|", lapply(as.list(...), function(n) grepl(glob2rx(n), names(x))))
  #   select.row <- Reduce("|", lapply(as.list(...), function(n) grepl(glob2rx(n), rownames(attr(x, "deriv")))))
  #   out <- unclass(x)[select.name]
  #   deriv <- submatrix(attr(x, "deriv"), rows = select.row)
  # } else {
  #   out <- unclass(x)[...]
  #   deriv <- submatrix(attr(x, "deriv"), rows = ...)
  # }
  # 
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


