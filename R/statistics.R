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
  sanePars <- sanitizePars(pars, list(...)$fixed)
  pars <- sanePars$pars
  fixed <- sanePars$fixed
  
  
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




#' Non-Linear Optimization, multi start
#' 
#' @description Wrapper around \code{\link{trust}} allowing for multiple fits 
#'   from randomly chosen initial values.
#'   
#' @param objfun Objective function, see \code{\link{trust}}.
#' @param center Parameter values around which the initial values for each fit 
#'   are randomly sampled. The initial values handed to \link{trust} are the sum
#'   of center and the output of \option{samplefun}, center + 
#'   \option{samplefun}. See \code{\link{trust}}, parinit.
#'   \code{center} Can also be a parframe, then the parameter values are taken 
#'   from the parframe. In this case, the \code{fits} argument is overwritten.
#' @param studyname The names of the study or fit. This name is used to 
#'   determine filenames for interim and final results. See Details.
#' @param rinit Starting trust region radius, see \code{\link{trust}}.
#' @param rmax Maximum allowed trust region radius, see \code{\link{trust}}.
#' @param fits Number of fits (jobs).
#' @param cores Number of cores for job parallelization.
#' @param samplefun Function to sample random initial values. It is assumed, 
#'   that \option{samplefun} has a named parameter "n" which defines how many 
#'   random numbers are to be returned, such as for \code{\link{rnorm}} or 
#'   \code{\link{runif}}. By default \code{\link{rnorm}} is used. Parameteres 
#'   for samplefun are simply appended as named parameters to the mstrust call 
#'   and automatically handed to samplefun by matching parameter names.
#' @param resultPath character indicating the folder where the results should 
#'   be stored. Defaults to ".". 
#' @param stats If true, the same summary statistic as written to the logfile is
#'   printed to command line on mstrust completion.
#' @param ... Additional parameters handed to trust(), samplefun(), or the 
#'   objective function by matching parameter names. All unmatched parameters 
#'   are handed to the objective function objfun(). The log file starts with a 
#'   table telling which parameter was assigend to which function.
#' @param output logical. If true, writes output to the disc.
#'   
#' @details By running multiple fits starting at randomly chosen inital 
#'   parameters, the chisquare landscape can be explored using a deterministic 
#'   optimizer. Here, \code{\link{trust}} is used for optimization. The standard
#'   procedure to obtain random initial values is to sample random variables 
#'   from a uniform distribution (\code{\link{rnorm}}) and adding these to 
#'   \option{center}. It is, however, possible, to employ any other sampling 
#'   strategy by handing the respective function to mstrust(), 
#'   \option{samplefun}.
#'   
#'   In case a special sampling is required, a customized sampling function can 
#'   be used. If, e.g., inital values leading to a non-physical systems are to 
#'   be discarded upfront, the objective function can be addapted accordingly.
#'   
#'   All started fits either lead to an error or complete converged or
#'   unconverged. A statistics about the return status of fits can be shown by
#'   setting \option{stats} to TRUE.
#'   
#'   Fit final and intermediat results are stored under \option{studyname}. For
#'   each run of mstrust for the same study name, a folder under
#'   \option{studyname} of the form "trial-x-date" is created. "x" is the number
#'   of the trial, date is the current time stamp. In this folder, the
#'   intermediate results are stored. These intermediate results can be loaded
#'   by \code{\link{load.parlist}}. These are removed on successfull completion
#'   of mstrust. In this case, the final list of fit parameters
#'   (parameterList.Rda) and the fit log (mstrust.log) are found instead.
#'   
#' @return A parlist holding errored and converged fits.
#'   
#' @seealso \code{\link{trust}}, \code{\link{rnorm}}, \code{\link{runif}}, 
#'   \code{\link{as.parframe}}
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'  
#' @export
#' @import parallel
mstrust <- function(objfun, center, studyname, rinit = .1, rmax = 10, fits = 20, cores = 1,
                    samplefun = "rnorm", resultPath = ".", stats = FALSE, output = FALSE,
                    ...) {

  narrowing <- NULL
  
  # Check if on Windows
  cores <- sanitizeCores(cores)
  
  # Argument parsing, sorting, and enhancing
  # Gather all function arguments
  varargslist <- list(...)

  argslist <- as.list(formals())
  argslist <- argslist[names(argslist) != "..."]
  
  argsmatch <- as.list(match.call(expand.dots = TRUE))
  namesinter <- intersect(names(argslist), names(argsmatch))

  argslist[namesinter] <- argsmatch[namesinter]
  argslist <- c(argslist, varargslist)
  
  argslist[["objfun"]] <- force(objfun)
  argslist[["center"]] <- force(center)

  # Add extra arguments
  argslist$n <- length(center) # How many inital values do we need?

  # Determine target function for each function argument.
  # First, define argument names used locally in mstrust().
  # Second, check what trust() and samplefun() accept and check for name clashes.
  # Third, whatever is unused is passed to the objective function objfun().
  nameslocal <- c("studyname", "center", "fits", "cores", "samplefun",
                  "resultPath", "stats", "narrowing")
  namestrust <- intersect(names(formals(trust)), names(argslist))
  namessample <- intersect(names(formals(samplefun)), names(argslist))
  if (length(intersect(namestrust, namessample) != 0)) {
    stop("Argument names of trust() and ", samplefun, "() clash.")
  }
  namesobj <- setdiff(names(argslist), c(namestrust, namessample, nameslocal))


  # Assemble argument lists common to all calls in mclapply
  # Sample function
  argssample <- structure(vector("list", length = length(namessample)), names = namessample)
  for (name in namessample) {
    argssample[[name]] <- argslist[[name]]
  }

  # Objective function
  argsobj <- structure(vector("list", length = length(namesobj)), names = namesobj)
  for (name in namesobj) {
    argsobj[[name]] <- argslist[[name]]
  }

  # Trust optimizer, except for initial values
  argstrust <- structure(vector("list", length = length(namestrust)), names = namestrust)
  for (name in namestrust) {
    argstrust[[name]] <- argslist[[name]]
  }


  # Assemble and create output filenames, folders and files
  m_timeStamp <- paste0(format(Sys.time(), "%d-%m-%Y-%H%M%S"))
  
  # Folders
  resultFolderBase <- file.path(argslist$resultPath, argslist$studyname)
  m_trial <- paste0("trial-", length(dir(resultFolderBase, pattern = "trial*")) + 1)
  resultFolder <- file.path(resultFolderBase, paste0(m_trial, "-", m_timeStamp))
  
  interResultFolder <- file.path(resultFolder, "interRes")
  dir.create(path = interResultFolder, showWarnings = FALSE, recursive = TRUE)
  
  # Files
  fileNameLog <- paste0("mstrust.log")
  fileNameParList <- paste0("parameterList.Rda")
  fileLog <- file.path(resultFolder, fileNameLog)
  fileParList <- file.path(resultFolder, fileNameParList)
  
  
  
  # Apply trust optimizer in parallel
  # The error checking leverages that mclappy runs each job in a try().
  logfile <- file(fileLog, open = "w")
  
  # Parameter assignment information
  if (is.null(narrowing) || narrowing[1] == 1) {
    msg <- paste0("Parameter assignment information\n",
                  strpad("mstrust", 12),                        ": ", paste0(nameslocal, collapse = ", "), "\n",
                  strpad("trust", 12),                          ": ", paste0(namestrust, collapse = ", "), "\n",
                  strpad(as.character(argslist$samplefun), 12), ": ", paste0(namessample, collapse = ", "), "\n\n")
                  #strpad(as.character(argslist$objfun), 12),    ": ", paste0(namesobj, collapse = ", "), "\n\n")
    writeLines(msg, logfile)
    flush(logfile)
  }

  # Write narrowing status information to file
  if (!is.null(narrowing)) {
    msg <- paste0("--> Narrowing, run ", narrowing[1], " of ", narrowing[2], "\n",
                  "--> " , fits, " fits to run\n")
    writeLines(msg, logfile)
    flush(logfile)
  }

  if(is.parframe(center)) {
    fits <- nrow(center)
  }
  
  m_parlist <- as.parlist(mclapply(1:fits, function(i) {
    
    if(is.parframe(center)) {
      argstrust$parinit <- as.parvec(center, i)
    } else {
      argstrust$parinit <- center + do.call(samplefun, argssample)
    }
    
    fit <- do.call(trust, c(argstrust, argsobj))

    # Keep only numeric attributes of object returned by trust()
    attr.fit <- attributes(fit)
    keep.attr <- sapply(attr.fit, is.numeric)
    fit <- fit[1:length(fit)] # deletes attributes
    if (any(keep.attr)) attributes(fit) <- c(attributes(fit), attr.fit[keep.attr]) # attach numeric attributes
    
    
    # In some crashes a try-error object is returned which is not a list. Since
    # each element in the parlist is assumed to be a list, we wrap these cases.
    if (!is.list(fit)) {
      f <- list()
      f$error <- fit
      fit <- f
    }
    
    fit$parinit <- argstrust$parinit

    # Write current fit to disk
    if (output) {
      saveRDS(fit, file = file.path(interResultFolder, paste0("fit-", i, ".Rda")))
      
      # Reporting
      # With concurent jobs and everyone reporting, this is a classic race
      # condition. Assembling the message beforhand lowers the risk of interleaved
      # output to the log.
      msgSep <- "-------"
      if (any(names(fit) == "error")) {
        msg <- paste0(msgSep, "\n",
                      "Fit ", i, " failed after ", fit$iterations, " iterations with error\n",
                      "--> ", fit$error,
                      msgSep, "\n")
        
        writeLines(msg, logfile)
        flush(logfile)
      } else {
        msg <- paste0(msgSep, "\n",
                      "Fit ", i, " completed\n",
                      "--> iterations : ", fit$iterations, "\n",
                      "-->  converged : ", fit$converged, "\n",
                      "--> obj. value : ", round(fit$value, digits = 2), "\n",
                      msgSep)
        
        writeLines(msg, logfile)
        flush(logfile)
      }
    }
    return(fit)
  }, mc.preschedule = FALSE, mc.silent = FALSE, mc.cores = cores))
  close(logfile)


  # Cull failed and completed fits Two kinds of errors occure. The first returns
  # an object of class "try-error". The reason for these failures are unknown to
  # me. The second returns a list of results from trust(), where one name of the
  # list is error holding an object of class "try-error". These abortions are 
  # due to errors which are captured within trust(). Completed fits return with 
  # a valid result list from trust(), with "error" not part of its names. These
  # fits, can still be unconverged, if the maximim number of iterations was the
  # reason for the return of trust(). Be also aware of fits which converge due
  # to the trust radius hitting rmin. Such fits are reported as converged but
  # are not in truth.
  m_trustFlags.converged = 0
  m_trustFlags.unconverged = 1
  m_trustFlags.error = 2
  m_trustFlags.fatal = 3
  idxStatus <- sapply(m_parlist, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error")) {
      return(m_trustFlags.error)
    } else if (!any(names(fit) == "converged")) {
      return(m_trustFlags.fatal)
    } else if (fit$converged) {
      return(m_trustFlags.converged)
    } else {
      return(m_trustFlags.unconverged)
    }
  })

  
  # Wrap up
  # Write out results
  if (output) saveRDS(m_parlist, file = fileParList)

  # Remove temporary files
  unlink(interResultFolder, recursive = TRUE)
  

  # Show summary
  sum.error <- sum(idxStatus == m_trustFlags.error)
  sum.fatal <- sum(idxStatus == m_trustFlags.fatal)
  sum.unconverged <- sum(idxStatus == m_trustFlags.unconverged)
  sum.converged <- sum(idxStatus == m_trustFlags.converged)
  msg <- paste0("Multi start trust summary\n",
                "Outcome     : Occurrence\n",
                "Error       : ", sum.error, "\n",
                "Fatal       : ", sum.fatal, " must be 0\n",
                "Unconverged : ", sum.unconverged, "\n",
                "Converged   : ", sum.converged, "\n",
                "           -----------\n",
                "Total       : ", sum.error + sum.fatal + sum.unconverged + sum.converged, paste0("[", fits, "]"), "\n")
  logfile <- file(fileLog, open = "a")
  writeLines(msg, logfile)
  flush(logfile)
  close(logfile)

  if (stats) {
    cat(msg)
  }

  return(m_parlist)
}



#' Construct fitlist from temporary files.
#'
#' @description An aborted \code{\link{mstrust}}
#'   leaves behind results of already completed fits. This command loads these
#'   fits into a fitlist.
#'
#' @param folder Path to the folder where the fit has left its results.
#'
#' @details The command \code{\link{mstrust}} saves
#'   each completed fit along the multi-start sequence such that the results can
#'   be resurected on abortion. This command loads a fitlist from these
#'   intermediate results.
#'
#' @return An object of class parlist.
#'
#' @seealso \code{\link{mstrust}}
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
load.parlist <- function(folder) {
  # Read in all fits
  m_fileList <- dir(folder, pattern = "*.Rda")
  m_parVec <- lapply(m_fileList, function(file) {
    return(readRDS(file.path(folder, file)))
  })
  
  return(as.parlist(m_parVec))
}




#' Reduce replicated measurements to mean and standard deviation
#'
#' @description
#' Obtain the mean and standard deviation from replicates per condition.
#'
#' @param file Data file of csv. See Format for details.
#' @param select Names of the columns in the data file used to define
#'        conditions, see Details.
#' @param datatrans Character vector describing a function to transform data.
#'        Use \kbd{x} to refere to data.
#'
#'
#' @format
#' The following columns are mandatory for the data file.
#' \describe{
#'  \item{name}{Name of the observed species.}
#'  \item{time}{Measurement time point.}
#'  \item{value}{Measurement value.}
#'  \item{condition}{The condition under which the observation was made.}
#' }
#'
#' In addition to these columns, any number of columns can follow to allow a
#' fine grained definition of conditions. The values of all columns named in
#' \option{select} are then merged to get the set of conditions.
#'
#' @details
#' Experiments are usually repeated multiple times possibly under different
#' conditions leading to replicted measurements. The column "Condition" in the
#' data allows to group the data by their condition. However, sometimes, a more
#' fine grained grouping is desirable. In this case, any number of additional
#' columns can be append to the data. These columns are referred to as
#' "condition identifier". Which of the condition identifiers are used to do the
#' grouping is user defined by anouncing the to \option{select}. The mandatory
#' column "Condition" is always used. The total set of different conditions is
#' thus defined by all combinations of values occuring in the selected condition
#' identifiers. The replicates of each condition are then reduced to mean and
#' variance.New conditions names are derived by merging all conditions which
#' were used in mean and std.
#'
#' @return
#' A data frame of the following variables
#' \describe{
#'  \item{time}{Measurement time point.}
#'  \item{name}{Name of the observed species.}
#'  \item{value}{Mean of replicates.}
#'  \item{sigma}{Standard error of the mean, NA for single measurements.}
#'  \item{n}{The number of replicates reduced.}
#'  \item{condition}{The condition for which the value and sigma were calculated. If
#'        more than one column were used to define the condition, this variable
#'        holds the effecive condition which is the combination of all applied
#'        single conditions. }
#' }
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
reduceReplicates <- function(file, select = "condition", datatrans = NULL) {

  # File format definiton
  fmtnames <- c("name", "time",  "value", "condition")
  fmtnamesnumber <- length(fmtnames)

  # Read data and sanity checks
  data <- read.csv(file)
  if (length(intersect(names(data), fmtnames)) != fmtnamesnumber) {
    stop(paste("Mandatory column names are:", paste(fmtnames, collapse = ", ")))
  }

  # Transform data if requested
  if (is.character(datatrans)) {
    x <- data$value
    data$value <- eval(parse(text = datatrans))
  }

  # Experiments are usually repeated multiple times possibly under different
  # conditions. The column "Condition" in the data thus groups the data per
  # condition. However, sometimes, a more fine grained grouping is desirable. In
  # this case, any number of additional columns can be append to the data. These
  # columns are referred to as "condition identifier". Which of the condition
  # identifiers are used to do the grouping is user defined by giving their
  # names in <select>. The mandatory column "Condition" is always used. The
  # total set of different conditions is thus defined by all combinations of
  # values occuring in the condition identifiers named for grouping. Mean and
  # variance is computed for each condition by averaging over measurements
  # recorded at the same time point. New conditions names are derived by merging
  # all conditions which were used in mean and std.
  select <- unique(c("name", "time", "condition", select))
  condidnt <- Reduce(paste, subset(data, select = select))
  conditions <- unique(condidnt)
  reduct <- do.call(rbind, lapply(conditions, function(cond) {
    conddata <- data[condidnt == cond,]
    mergecond <- Reduce(paste, conddata[1, setdiff(select, c("name", "time"))])
    data.frame(time = conddata[1, "time"],
               value = mean(conddata[, "value"]),
               sigma = sd(conddata[, "value"])/sqrt(nrow(conddata)),
               n = nrow(conddata),
               name = conddata[1, "name"],
               condition = mergecond)
  }))

  return(reduct)
}



#' Fit an error model
#'
#' @description Fit an error model to reduced replicate data, see
#'   \code{\link{reduceReplicates}}.
#'
#' @param data Reduced replicate data, see \code{\link{reduceReplicates}}. Need 
#'   columns "value", "sigma", "n".
#' @param factors \option{data} is pooled with respect to the columns named
#'   here, see Details.
#' @param errorModel Character vector defining the error model in terms of the variance. 
#'   Use \kbd{x} to reference the independend variable, see Details.
#' @param par Inital values for the parameters of the error model.
#' @param plotting If TRUE, a plot of the pooled variance together with the fit
#'   of the error model is shown.
#' @param blather If TRUE, additional information is returned, such as fit parameters 
#'  and sigmaLS (original sigma given in input data).
#' @param ... Parameters handed to the optimizer \code{\link{optim}}.
#'
#' @details The variance estimator using \eqn{n-1} data points is \eqn{chi^2}
#'   distributed with \eqn{n-1} degrees of freedom. Given replicates for
#'   consecutive time points, the sample variance can be assumed a function of
#'   the sample mean. By defining an error model which must hold for all time
#'   points, a maximum likelihood estimator for the parameters of the error
#'   model can be derived. The parameter \option{errorModel} takes the error
#'   model as a character vector, where the mean (independent variable) is
#'   refered to as \kbd{x}.
#'
#'   It is desireable to estimate the variance from many replicates. The
#'   parameter \option{data} must provide one or more columns which define the
#'   pooling of data. In case more than one column is announced by
#'   \option{factors}, all combinations are constructed. If, e.g.,
#'   \option{factors = c("condition", "name")} is used, where "condition" is
#'   "a", "b", "c" and repeating and "name" is "d", "e" and repeating, the
#'   effective conditions used for pooling are "a d", "b e", "c d", "a e", "b
#'   d", and "c e".
#'
#'   By default, a plot of the pooled data, sigma and its confidence bound at
#'   68\% and 95\% is shown.
#'
#' @return Returned by default is a data frame with columns as in \option{data}, 
#'   but with the sigma values replaced by the derived values, obtained by evaluating 
#'   the error model with the fit parameters. 
#'   
#'   If the blather = TRUE option is chosen, fit values of the parameters of the error
#'   model are appended, with the column names equal to the parameter names. 
#'   The error model is appended as the attribute "errorModel".
#'   Confidence bounds for sigma at confidence level 68\% and 95\% are
#'   calculated, their values come next in the returned data frame. Finally, the
#'   effective conditions are appended to easily check how the pooling was done.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
#' @importFrom stats D approx optim qchisq sd time
fitErrorModel <- function(data, factors, errorModel = "exp(s0)+exp(srel)*x^2",
                          par = c(s0 = 1, srel = .1), plotting = TRUE, blather = FALSE, ...) {

  # Assemble conditions
  condidnt <- Reduce(paste, subset(data, select = factors))
  conditions <- unique(condidnt)


  # Fit error model
  nColData <- ncol(data)
  dataErrorModel <- cbind(data, as.list(par))

  for (cond in conditions) {
    subdata <- dataErrorModel[condidnt == cond,]
    x <- subdata$value
    n <- subdata$n
    y <- subdata$sigma*sqrt(n)

    obj <- function(par) {
      value <- with(as.list(par), {
        z <- eval(parse(text = errorModel))
        sum(log(z)-log(dchisq((n-1)*(y^2)/z, df = n-1)), na.rm = TRUE)
      })
      return(value)
    }

    fit <- optim(par = par, fn = obj, ...)
    sigma <- sqrt(with(as.list(fit$par), eval(parse(text = errorModel))))
    dataErrorModel[condidnt == cond, ]$sigma <- sigma 
    dataErrorModel[condidnt == cond, -(nColData:1)] <- data.frame(as.list(fit$par))
  }


  # Calculate confidence bounds about sigma
  p68 <- (1-.683)/2
  p95 <- (1-.955)/2
  dataErrorModel$cbLower68 <- dataErrorModel$sigma^2*qchisq(p = p68, df = dataErrorModel$n-1)/(dataErrorModel$n-1)
  dataErrorModel$cbUpper68 <- dataErrorModel$sigma^2*qchisq(p = p68, df = dataErrorModel$n-1, lower.tail = FALSE)/(dataErrorModel$n-1)
  dataErrorModel$cbLower95 <- dataErrorModel$sigma^2*qchisq(p = p95, df = dataErrorModel$n-1)/(dataErrorModel$n-1)
  dataErrorModel$cbUpper95 <- dataErrorModel$sigma^2*qchisq(p = p95, df = dataErrorModel$n-1, lower.tail = FALSE)/(dataErrorModel$n-1)


  # Assemble result
  dataErrorModel <- cbind(dataErrorModel, condidnt, sigmaLS = data$sigma)
  attr(dataErrorModel, "errorModel") <- errorModel


  # Plot if requested
  if (plotting) {
    print(ggplot(dataErrorModel, aes(x=value)) +
            geom_point(aes(y=sigmaLS^2*(n))) +
            geom_line(aes(y=sigma^2)) +
            geom_ribbon(aes(ymin=cbLower95, ymax=cbUpper95), alpha=.3) +
            geom_ribbon(aes(ymin=cbLower68, ymax=cbUpper68), alpha=.3) +
            ylab("variance") +
            facet_wrap(~condidnt, scales = "free") +
            scale_y_log10() +
            theme_dMod()
    )}
  
  # Return standard error of the mean
  dataErrorModel$sigma <- dataErrorModel$sigma/sqrt(dataErrorModel$n)
  data$sigma <- dataErrorModel$sigma
  if(blather)
    return(dataErrorModel)
  else 
    return(data)
}




