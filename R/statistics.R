#' Profile-likelihood (PL) computation
#' 
#' @param obj Objective function \code{obj(pars, fixed, ...)} returning a list with "value",
#' "gradient" and "hessian". If attribute "valueData" and/or "valuePrior are returned they are attached to the return value.
#' @param pars Parameter vector corresponding to the log-liklihood optimum.
#' @param whichPar Numeric or character vector. The parameters for which the profile is computed.
#' @param alpha Numeric, the significance level based on the chisquare distribution with df=1
#' @param limits Numeric vector of length 2, the lower and upper deviance from the original 
#' value of \code{pars[whichPar]}
#' @param method Character, either \code{"integrate"} or \code{"optimize"}. This is a short-cut for
#' setting stepControl, algoControl and optControl by hand.
#' @param stepControl List of arguments controlling the step adaption. Defaults to integration set-up, i.e.
#' \code{list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 100)}
#' @param algoControl List of arguments controlling the fast PL algorithm. defaults to
#' \code{list(gamma = 1, W = "hessian", reoptimize = FALSE, correction = 1, reg = .Machine$double.eps)}
#' @param optControl List of arguments controlling the \code{trust()} optimizer. Defaults to
#' \code{list(rinit = .1, rmax = 10, iterlim = 10, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))}.
#' See \link{trust} for more details.
#' @param verbose Logical, print verbose messages.
#' @param cores number of cores used when computing profiles for several
#' parameters.
#' @param ... Arguments going to obj()
#' @details Computation of the profile likelihood is based on the method of Lagrangian multipliers
#' and Euler integration of the corresponding differential equation of the profile likelihood paths.
#' 
#' \code{algoControl}: Since the Hessian which is needed for the differential equation is frequently misspecified, 
#' the error in integration needs to be compensated by a correction factor \code{gamma}. Instead of the
#' Hessian, an identity matrix can be used. To guarantee that the profile likelihood path stays on
#' the true path, each point proposed by the differential equation can be used as starting point for
#' an optimization run when \code{reoptimize = TRUE}. The correction factor \code{gamma} is adapted
#' based on the amount of actual correction. If this exceeds the value \code{correction}, \code{gamma} is
#' reduced. In some cases, the Hessian becomes singular. This leads to problems when inverting the
#' Hessian. To avoid this problem, the pseudoinverse is computed by removing all singular values lower
#' than \code{reg}.
#' 
#' \code{stepControl}: The Euler integration starts with \code{stepsize}. In each step the predicted change
#' of the objective function is compared with the actual change. If this is larger than \code{atol}, the
#' stepsize is reduced. For small deviations, either compared the abolute tolerance \code{atol} or the
#' relative tolerance \code{rtol}, the stepsize may be increased. \code{max} and \code{min} are upper and lower
#' bounds for \code{stepsize}. \code{limit} is the maximum number of steps that are take for the profile computation.
#' \code{stop} is a character, usually "value" or "data", for which the significance level \code{alpha}
#' is evaluated.
#' 
#' @return Named list of length one. The name is the parameter name. The list enty is a
#' matrix with columns "value" (the objective value), "constraint" (deviation of the profiled paramter from
#' the original value), "stepsize" (the stepsize take for the iteration), "gamma" (the gamma value employed for the
#' iteration), "valueData" and "valuePrior" (if specified in obj), one column per parameter (the profile paths).
#' @example inst/examples/profiles.R
#' @export
profile <- function(obj, pars, whichPar, alpha = 0.05, 
                          limits = c(lower = -Inf, upper = Inf), 
                          method = c("integrate", "optimize"),
                          stepControl = NULL, 
                          algoControl = NULL,
                          optControl  = NULL,
                          verbose = FALSE,
                          cores = 1,
                          ...) {
  
  # Guarantee that pars is named numeric without deriv attribute
  pars <- structure(as.numeric(pars), names = names(pars))
  
  # Initialize control parameters depending on method
  method  <- match.arg(method)
  if (method == "integrate") {
    sControl <- list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 500, stop = "value")
    aControl <- list(gamma = 1, W = "hessian", reoptimize = FALSE, correction = 1, reg = .Machine$double.eps)
    oControl <- list(rinit = .1, rmax = 10, iterlim = 10, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  if (method == "optimize") {
    sControl <- list(stepsize = 1e-2, min = 1e-4, max = Inf, atol = 1e-1, rtol = 1e-1, limit = 100, stop = "value")
    aControl <- list(gamma = 0, W = "identity", reoptimize = TRUE, correction = 1, reg = 0)
    oControl <- list(rinit = .1, rmax = 10, iterlim = 100, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  
  # Check if on Windows
  cores <- sanitizeCores(cores)
  
  # Substitute user-set control parameters
  if (!is.null(stepControl)) sControl[match(names(stepControl), names(sControl))] <- stepControl
  if (!is.null(algoControl)) aControl[match(names(algoControl), names(aControl))] <- algoControl
  if (!is.null(optControl )) oControl[match(names(optControl), names(oControl ))] <- optControl
    
  do.call(rbind, mclapply(whichPar, function(whichPar) {
    
    
    if (is.character(whichPar)) whichPar <- which(names(pars) == whichPar)
    whichPar.name <- names(pars)[whichPar]
    
    
    
    
    if (any(names(list(...)) == "fixed")) fixed <- list(...)$fixed else fixed <- NULL
    
    
    ## Functions needed during profile computation -----------------------
    obj.opt <- obj
    obj.prof <- function(p, ...) {
      out <- obj(p, ...)
      # If "identity", substitute hessian such that steps are in whichPar-direction.
      Id <- diag(1/.Machine$double.eps, length(out$gradient))
      Id[whichPar, whichPar] <- 1
      colnames(Id) <- rownames(Id) <- names(out$gradient)
      
      W <- match.arg(aControl$W[1], c("hessian", "identity"))
      out$hessian <- switch(W,
                            "hessian" = out$hessian,
                            "identity" = Id)
      return(out)    
    }
    
    pseudoinverse <- function(m, tol) {
      msvd <- svd(m)
      index <- which(abs(msvd$d) > max(dim(m))*max(msvd$d)*tol) 
      if (length(index) == 0) {
        out <- array(0, dim(m)[2:1])
      }
      else {
        out <- msvd$u[,index] %*% (1/msvd$d[index] * t(msvd$v)[index,])
      }
      attr(out, "valid") <- 1:length(msvd$d) %in% index
      return(out)
    }
    
    constraint <- function(p) {
      value <- p[whichPar] - pars[whichPar]
      gradient <- rep(0, length(p))
      gradient[whichPar] <- 1
      return(list(value = value, gradient = gradient))
    }
    lagrange <- function(y) {
      
      # initialize values
      p <- y
      lambda <- 0
      out <- obj.prof(p, ...)
      g.original <- constraint(p)
      
      # evaluate derivatives and constraints
      g     <- direction * g.original$value
      gdot  <- direction * g.original$gradient
      ldot  <- out$gradient
      lddot <- out$hessian 
      
      # compute rhs of profile ODE
      M <- rbind(cbind(lddot, gdot), 
                 matrix(c(gdot, 0), nrow=1))
      
      v <- c(-rep(gamma, length(p))*(ldot + lambda*gdot), 1)
      v0 <- c(-rep(0, length(p))*(ldot + lambda*gdot), 1)
      
      W <- pseudoinverse(M, tol = aControl$reg)
      valid <- attr(W, "valid")
      if(any(!valid)) {
        dy <- try(as.vector(W%*%v)[1:length(p)], silent=FALSE)
        dy0 <- try(as.vector(W%*%v0)[1:length(p)], silent=FALSE)
        dy[!valid[1:length(p)]] <- dy0[!valid[1:length(p)]] <- 0
        dy[whichPar] <- dy0[whichPar] <- direction
        warning(paste0("Iteration ", i, ": Some singular values of the Hessian are below the threshold. Optimization will be performed."))
      } else {
        dy <- try(as.vector(W%*%v)[1:length(p)], silent=FALSE)
        dy0 <- try(as.vector(W%*%v0)[1:length(p)], silent=FALSE)
      }
      
      
      if(!inherits(dy, "try-error")) {
        names(dy) <- names(y) 
        correction <- sqrt(sum((dy-dy0)^2))/sqrt(sum(dy^2))
      } else {
        dy <- NA
        correction <- 0
        warning(paste0("Iteration ", i, ": Impossible to invert Hessian. Trying to optimize instead."))
      }
      
      # Get attributes of the returned objective value (only if numeric)
      out.attributes <- attributes(out)[sapply(attributes(out), is.numeric)]
      out.attributes.names <- names(out.attributes)
      
      
      return(c(list(dy = dy, value = out$value, gradient = out$gradient, correction = correction, valid = valid, attributes = out.attributes.names),
               out.attributes))
      
      
      
      
    }
    doIteration <- function() {
      
      optimize <- aControl$reoptimize
      # Check for error in evaluation of lagrange()
      if(is.na(dy[1])) {
        #cat("Evaluation of lagrange() not successful. Will optimize instead.\n")
        optimize <- TRUE
        y.try <- y
        y.try[whichPar] <- y[whichPar] + direction*stepsize
        rinit <- oControl$rinit
      } else {
        dy.norm <- sqrt(sum(dy^2))
        rinit <- min(c(oControl$rinit, 3*dy.norm))
        y.try <- y + dy
        if(any(!lagrange.out$valid)) optimize <- TRUE
      }
      
      # Do reoptimization if requested or necessary
      if(optimize) {      
        parinit.opt <- y.try[-whichPar]
        fixed.opt <- c(fixed, y.try[whichPar])
        
        arglist <- c(list(objfun = obj.opt, parinit = parinit.opt, fixed = fixed.opt, rinit = rinit), 
                     oControl[names(oControl)!="rinit"],
                     list(...)[names(list(...)) != "fixed"])
        
        
        myfit <- try(do.call(trust, arglist), silent=FALSE)
        if(!inherits(myfit, "try-error")) {
          y.try[names(myfit$argument)] <- as.vector(myfit$argument)  
        } else {
          warning("Optimization not successful. Profile may be erroneous.")
        }
        
      }
      
      return(y.try)
      
    }
    doAdaption <- function() {
      
      lagrange.out.try <- lagrange(y.try)
      valid <- TRUE
      
      # Predicted change of the objective value
      dobj.pred <- sum(lagrange.out$gradient*(y.try - y))
      dobj.fact <- lagrange.out.try$value - lagrange.out$value
      correction <- lagrange.out.try$correction
      
      # Gamma adaption based on amount of actual correction
      if (correction > aControl$correction) gamma <- gamma/2
      if (correction < 0.5*aControl$correction) gamma <- min(c(aControl$gamma, gamma*2))
      
      # Stepsize adaption based on difference in predicted change of objective value
      if (abs(dobj.fact - dobj.pred) > sControl$atol & stepsize > sControl$min) {
        stepsize <- max(c(stepsize/1.5, sControl$min))
        valid <- FALSE
      }
      if (abs(dobj.fact - dobj.pred) < .3*sControl$atol | abs((dobj.fact - dobj.pred)/dobj.fact) < .3*sControl$rtol) {
        stepsize <- min(c(stepsize*2, sControl$max))
      }
      
      ## Verbose
      if (verbose) {
        # Compute progres
        diff.thres <- diff.steps <- diff.limit <- 0
        if (threshold < Inf)
          diff.thres <- 1 - max(c(0, min(c(1, (threshold - lagrange.out.try$value)/delta))))
        if (sControl$limit < Inf)
          diff.steps <- i/sControl$limit
        diff.limit <- switch(as.character(sign(constraint.out$value)),
                             "1"  = 1 - (limits[2] - constraint.out$value)/limits[2],
                             "-1" = diff.limit <- 1 - (limits[1] - constraint.out$value)/limits[1],
                             "0"  = 0)
        
        percentage <- max(c(diff.thres, diff.steps, diff.limit), na.rm = TRUE)*100
        progressBar(percentage)
        
        
        #cat("diff.thres:", diff.thres, "diff.steps:", diff.steps, "diff.limit:", diff.limit)
        myvalue <- format(substr(lagrange.out$value  , 0, 8), width = 8)
        myconst <- format(substr(constraint.out$value, 0, 8), width = 8)
        mygamma <- format(substr(gamma               , 0, 8), width = 8)
        myvalid <- all(lagrange.out$valid)
        cat("\tvalue:", myvalue, "constraint:", myconst, "gamma:", mygamma, "valid:", myvalid) 
      }
      
      
      
      return(list(lagrange = lagrange.out.try, stepsize = stepsize, gamma = gamma, valid = valid))
      
      
    }
    
    ## Compute profile -------------------------------------------------
    
    # Initialize profile
    i <- 0 
    direction <- 1
    gamma <- aControl$gamma
    stepsize <- sControl$stepsize
    ini <- pars
    
    lagrange.out <- lagrange(ini)
    constraint.out <- constraint(pars)
    
    delta <- qchisq(1-alpha, 1)
    #delta.t <- qf(1 - alpha, 1L, nobs - npar)
    
    
    
    threshold <- lagrange.out[[sControl$stop]] + delta
    out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
    
    out <- c(value = lagrange.out$value, 
             constraint = as.vector(constraint.out$value), 
             stepsize = stepsize, 
             gamma = gamma, 
             whichPar = whichPar,
             out.attributes, ini)
    
    # Compute right profile
    if (verbose) {
      cat("Compute right profile\n")
    }
    direction <- 1
    gamma <- aControl$gamma
    stepsize <- sControl$stepsize
    y <- ini
    
    lagrange.out <- lagrange.out
    constraint.out <- constraint.out
    
    while (i < sControl$limit) {
      
      ## Iteration step
      sufficient <- FALSE
      retry <- 0
      while (!sufficient & retry < 5) {
        dy <- stepsize*lagrange.out$dy
        y.try <- try(doIteration(), silent = TRUE)
        out.try <- try(doAdaption(), silent = TRUE)
        if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) {
          sufficient <- FALSE
          stepsize <- stepsize/1.5
          retry <- retry + 1
        } else {
          sufficient <- out.try$valid
          stepsize <- out.try$stepsize  
        }
        
      }    
      if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) break
      
      
      ## Set values
      y <- y.try
      lagrange.out <- out.try$lagrange
      constraint.out <- constraint(y.try)
      stepsize <- out.try$stepsize
      gamma <- out.try$gamma
      out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
      
      ## Return values 
      out <- rbind(out, 
                   c(value = lagrange.out$value, 
                     constraint = as.vector(constraint.out$value), 
                     stepsize = stepsize, 
                     gamma = gamma, 
                     whichPar = whichPar,
                     out.attributes, 
                     y))
      
      value <- lagrange.out[[sControl$stop]]
      if (value > threshold | constraint.out$value > limits[2]) break
      
      i <- i + 1
      
    }
    
    # Compute left profile
    if (verbose) {
      cat("\nCompute left profile\n")
    }
    i <- 0
    direction <- -1
    gamma <- aControl$gamma
    stepsize <- sControl$stepsize
    y <- ini
    
    lagrange.out <- lagrange(ini)
    constraint.out <- constraint(pars)
    
    while (i < sControl$limit) {
      
      ## Iteration step
      sufficient <- FALSE
      retry <- 0
      while (!sufficient & retry < 5) {
        dy <- stepsize*lagrange.out$dy
        y.try <- try(doIteration(), silent = TRUE)
        out.try <- try(doAdaption(), silent = TRUE)
        if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) {
          sufficient <- FALSE
          stepsize <- stepsize/1.5
          retry <- retry + 1
        } else {
          sufficient <- out.try$valid
          stepsize <- out.try$stepsize  
        }
        
      }
      if (inherits(y.try, "try-error") | inherits(out.try, "try-error")) break
      
      ## Set values
      y <- y.try
      lagrange.out <- out.try$lagrange
      constraint.out <- constraint(y.try)
      stepsize <- out.try$stepsize
      gamma <- out.try$gamma
      out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
      
      
      ## Return values
      out <- rbind(c(value = lagrange.out$value, 
                     constraint = as.vector(constraint.out$value), 
                     stepsize = stepsize, 
                     gamma = gamma,
                     whichPar = whichPar,
                     out.attributes,
                     y), 
                   out)
      
      value <- lagrange.out[[sControl$stop]]
      if (value > threshold | constraint.out$value < limits[1]) break
      
      i <- i + 1
      
    }
    
    # Output
    out <- as.data.frame(out)
    out$whichPar <- whichPar.name
    parframe(
      out,
      parameters = names(pars),
      metanames = c("value", "constraint", "stepsize", "gamma", "whichPar"),
      obj.attributes = names(out.attributes)
    )
    
  }, mc.cores = cores, mc.preschedule = FALSE))  
  
}

#' Progress bar
#' 
#' @param percentage Numeric between 0 and 100
#' @param size Integer, the size of the bar print-out
#' @param number Logical, Indicates whether the percentage should be printed out.
#' @export
progressBar <- function(percentage, size = 50, number = TRUE) {
  
  if(percentage < 0) percentage <- 0
  if(percentage > 100) percentage <- 100
  
  out <- paste("\r|", paste(rep("=", round(size*percentage/100)), collapse=""), paste(rep(" ", size-round(size*percentage/100)), collapse=""), "|", sep="")
  cat(out)
  if(number) cat(format(paste(" ", round(percentage), "%", sep=""), width=5))
  
}

#' Profile uncertainty extraction
#' 
#' @description extract parameter uncertainties from profiles
#' @param object object of class \code{parframe}, returned from \link{profile} function.
#' @param parm a specification of which parameters are to be given confidence intervals, 
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... not used right now.
#' @param val.column the value column used in the parframe, usually 'data'.
#' @export
confint.parframe <- function(object, parm = NULL, level = 0.95, ..., val.column = "data") {
  
  profile <- object
  
  maxvalue <- qchisq(level, df = 1)
  
  proflist <- as.data.frame(profile)
  obj.attributes <- attr(profile, "obj.attributes")
  
  if(is.data.frame(proflist)) {
    whichPars <- unique(proflist$whichPar)
    if (!is.null(parm)) whichPars <- intersect(whichPars, parm)
    if (length(whichPars) == 0) stop("profile does not contain the required parameters.")
    proflist <- lapply(whichPars, function(n) {
      with(proflist, proflist[whichPar == n, ])
    })
    names(proflist) <- whichPars
  }
  
  # Discard faulty profiles
  proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
  proflist <- proflist[proflistidx]
  if (sum(!proflistidx) > 0) {
    warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
  }
  subdata <- do.call(rbind, lapply(names(proflist), function(n) {
    #print(n)
    parvalues <- proflist[[n]][, n]
    values <- proflist[[n]][, val.column]
    origin <- which.min(abs(proflist[[n]][, "constraint"]))
    zerovalue <- proflist[[n]][origin, val.column]
    deltavalues <- values - zerovalue
    deltavalues[deltavalues < 0] <- 0
    deltavalues <- sqrt(deltavalues)
    deltavalues[1:origin] <- - deltavalues[1:origin]
    lpars <- length(parvalues)
    x <- abs(deltavalues - sqrt(maxvalue))
    position_upper <- which(x %in% sort(x)[1:2])[1:2]
    x <- abs(deltavalues + sqrt(maxvalue))
    position_lower <- which(x %in% sort(x)[1:2])[1:2]
    #cat("deltas:",deltavalues[position_lower],deltavalues[position_upper], "\n")
    #cat("pars:",parvalues[position_lower],parvalues[position_upper], "\n")
    parlower <- try(approxExtrap(deltavalues[position_lower], y = parvalues[position_lower], xout = -sqrt(maxvalue)), silent = TRUE)
    parupper <- try(approxExtrap(deltavalues[position_upper], y = parvalues[position_upper], xout = sqrt(maxvalue)), silent = TRUE)
    
    if (inherits(parlower, "try-error")) parlower <- list(x = NA, y = NA)
    if (inherits(parupper, "try-error")) parupper <- list(x = NA, y = NA)
    
    parmin <- parvalues[origin]
    data.frame(name = n, value = parmin, lower = parlower$y, upper = parupper$y)
  }))
  return(subdata)
}

# Extrapolation extension to approx
# 
# from Hmisc package, built on approx
approxExtrap <- function (x, y, xout, method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE) {
if (is.list(x)) {
  y <- x[[2]]
  x <- x[[1]]
}
if (na.rm) {
  d <- !is.na(x + y)
  x <- x[d]
  y <- y[d]
}
d <- !duplicated(x)
x <- x[d]
y <- y[d]
d <- order(x)
x <- x[d]
y <- y[d]
w <- approx(x, y, xout = xout, method = method, n = n, rule = 2, 
            f = f, ties = ties)$y
r <- range(x)
d <- xout < r[1]
if (any(is.na(d))) 
  stop("NAs not allowed in xout")
if (any(d)) 
  w[d] <- (y[2] - y[1])/(x[2] - x[1]) * (xout[d] - x[1]) + 
  y[1]
d <- xout > r[2]
n <- length(y)
if (any(d)) 
  w[d] <- (y[n] - y[n - 1])/(x[n] - x[n - 1]) * (xout[d] - 
                                                   x[n - 1]) + y[n - 1]
list(x = xout, y = w)
}