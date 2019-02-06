
norm <- function(x) sqrt(sum(x^2))


#' Non-Linear Optimization
#' 
#' This function carries out a minimization or maximization of a function 
#' using a trust region algorithm. See the references for details.
#' 
#' @param objfun an R function that computes value, gradient, and Hessian of the 
#' function to be minimized or maximized and returns them as a list with 
#' components value, gradient, and hessian. Its first argument should be a 
#' vector of the length of parinit followed by any other arguments specified 
#' by the \code{...} argument.
#' 
#' @param parinit starting parameter values for the optimization. Must be 
#' feasible (in the domain).
#' 
#' @param rinit starting trust region radius. The trust region radius 
#' (see details below) is adjusted as the algorithm proceeds. A bad initial 
#' value wastes a few steps while the radius is adjusted, but does not keep 
#' the algorithm from working properly.
#' 
#' @param rmax maximum allowed trust region radius. This may be set very large. 
#' If set small, the algorithm traces a steepest descent path (steepest ascent, 
#' when minimize = FALSE).
#' 
#' @param parscale an estimate of the size of each parameter at the minimum. 
#' The algorithm operates as if optimizing function(x, ...) objfun(x / parscale, ...). 
#' May be missing in which case no rescaling is done. See also the details section below.
#' 
#' @param iterlim a positive integer specifying the maximum number of iterations 
#' to be performed before the program is terminated.
#' 
#' @param fterm a positive scalar giving the tolerance at which the difference 
#' in objective function values in a step is considered close enough to zero to 
#' terminate the algorithm.
#' 
#' @param mterm a positive scalar giving the tolerance at which the two-term 
#' Taylor-series approximation to the difference in objective function values in 
#' a step is considered close enough to zero to terminate the algorithm.
#' 
#' @param minimize If TRUE minimize. If FALSE maximize.
#' 
#' @param blather If TRUE return extra info.
#' 
#' @param parupper named numeric vector of upper bounds.
#' @param parlower named numeric vector of lower bounds.
#' 
#' @param printIter print iteration information to R console
#' 
#' @param ... additional argument to objfun
#' 
#' @details See Fletcher (1987, Section 5.1) or Nocedal and Wright (1999, Section 4.2) 
#' for detailed expositions.
#' 
#' @return A list containing the following components:
#' \itemize{
#' \item{value: }{the value returned by objfun at the final iterate.}
#' \item{gradient: }{the gradient returned by objfun at the final iterate.}
#' \item{hessian: }{the Hessian returned by objfun at the final iterate.}
#' \item{argument: }{the final iterate}
#' \item{converged: }{if TRUE the final iterate was deemed optimal by the 
#' specified termination criteria.}
#' \item{iterations: }{number of trust region subproblems done (including those 
#' whose solutions are not accepted).}
#' \item{argpath: }{(if blather == TRUE) the sequence of iterates, not including 
#' the final iterate.}
#' \item{argtry: }{(if blather == TRUE) the sequence of solutions of the trust 
#' region subproblem.}
#' \item{steptype: }{(if blather == TRUE) the sequence of cases that arise in 
#' solutions of the trust region subproblem. "Newton" means the Newton step 
#' solves the subproblem (lies within the trust region). Other values mean the 
#' subproblem solution is constrained. "easy-easy" means the eigenvectors 
#' corresponding to the minimal eigenvalue of the rescaled Hessian are not all 
#' orthogonal to the gradient. The other cases are rarely seen. "hard-hard" means 
#' the Lagrange multiplier for the trust region constraint is minus the minimal 
#' eigenvalue of the rescaled Hessian; "hard-easy" means it isn't.}
#' \item{accept: }{(if blather == TRUE) indicates which of the sequence of 
#' solutions of the trust region subproblem were accepted as the next iterate. 
#' (When not accepted the trust region radius is reduced, and the previous iterate 
#' is kept.)}
#' \item{r: }{(if blather == TRUE) the sequence of trust region radii.}
#' \item{rho: }{(if blather == TRUE) the sequence of ratios of actual over 
#' predicted decrease in the objective function in the trust region subproblem, 
#' where predicted means the predicted decrease in the two-term Taylor series model 
#' used in the subproblem.}
#' \item{valpath: }{(if blather == TRUE) the sequence of objective function values 
#' at the iterates.}
#' \item{valtry: }{(if blather == TRUE) the sequence of objective function values 
#' at the solutions of the trust region subproblem.}
#' \item{preddiff: }{(if blather == TRUE) the sequence of predicted differences using 
#' the two-term Taylor-series model between the function values at the current iterate 
#' and at the solution of the trust region subproblem.}
#' \item{stepnorm: }{(if blather == TRUE) the sequence of norms of steps, that is 
#' distance between current iterate and proposed new iterate found in the trust region 
#' subproblem.}
#' }
#' 
#' @export
#' @importFrom stats uniroot
trust <- function(objfun, parinit, rinit, rmax, parscale, iterlim = 100, 
                  fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps), 
                  minimize = TRUE, blather = FALSE, parupper = Inf, parlower = -Inf, printIter = FALSE, ...) 
{
  # Initialize ----
  # Guarantee that pars is named numeric without deriv attribute
  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit <- sanePars$pars
  
  
  u <- structure(rep(Inf, length(parinit)), names = names(parinit))
  l <- structure(rep(-Inf, length(parinit)), names = names(parinit))
  
  if (is.null(names(parupper)))
    u[1:length(u)] <- parupper
  if (is.null(names(parlower)))
    l[1:length(l)] <- parlower
  if (!is.null(names(parupper)) & !is.null(names(parinit)))
    u[names(parupper)] <- parupper
  if (!is.null(names(parlower)) & !is.null(names(parinit)))
    l[names(parlower)] <- parlower
  
  parupper <- u
  parlower <- l
  
  
  
  if (!is.numeric(parinit)) 
    stop("parinit not numeric")
  if (!all(is.finite(parinit))) 
    stop("parinit not all finite")
  
  upper <- which(parinit > parupper)
  lower <- which(parinit < parlower)
  if(length(upper > 0)){
    warning("init above range")
    parinit[upper] <- parupper[upper]
  }
  if(length(lower > 0)){
    warning("init below range")
    parinit[lower] <- parlower[lower]
  }
  d <- length(parinit)
  if (missing(parscale)) {
    rescale <- FALSE
  }
  else {
    rescale <- TRUE
    if (length(parscale) != d) 
      stop("parscale and parinit not same length")
    if (!all(parscale > 0)) 
      stop("parscale not all positive")
    if (!all(is.finite(parscale) & is.finite(1/parscale))) 
      stop("parscale or 1 / parscale not all finite")
  }
  if (!is.logical(minimize)) 
    stop("minimize not logical")
  r <- rinit
  theta <- parinit
  out <- try(objfun(theta, ...))
  if (inherits(out, "try-error")) {
    warning("error in first call to objfun")
    return(list(error = out, argument = theta, converged = FALSE, 
                iterations = 0))
  }
  checks <- try(check.objfun.output(out, minimize, d))
  if (inherits(checks, "try-error")) {
    warning("error in first call to objfun")
    return(c(as.list(out), error = checks, list(argument = theta), converged = FALSE, 
             iterations = 0))
  }
  if (!is.finite(out$value))  {
    error <- try(stop("parinit not feasible: value is not finite"))
    return(c(as.list(out), error = error, list(argument = theta), converged = FALSE, 
             iterations = 0))
  }
    
  
  #remove boundary elements from gradient and hessian
  # g_boundary <- c(upper[which(out$gradient[upper] < 0)], lower[which(out$gradient[lower] > 0)])
  # g_noboundary <- setdiff(1:d, g_boundary)
  # n_boundary <- length(g_boundary)
  # if (n_boundary > 0) {
  #   out$gradient <- out$gradient[-g_boundary]
  #   out$hessian <- out$hessian[-g_boundary, , drop = FALSE]
  #   out$hessian <- out$hessian[, -g_boundary, drop = FALSE]
  # }
  
  accept <- TRUE
  if (blather) {
    theta.blather <- NULL
    theta.try.blather <- NULL
    type.blather <- NULL
    accept.blather <- NULL
    r.blather <- NULL
    stepnorm.blather <- NULL
    rho.blather <- NULL
    val.blather <- NULL
    val.try.blather <- NULL
    preddiff.blather <- NULL
  }
  # Iterate ----
  
  if (printIter) cat("\n")
  
  for (iiter in 1:iterlim) {
    #cat(iiter, out$value,upper,"\n")
    if (blather) {
      theta.blather <- rbind(theta.blather, theta)
      r.blather <- c(r.blather, r)
      if (accept) 
        val.blather <- c(val.blather, out$value)
      else val.blather <- c(val.blather, out.value.save)
    }
    
    if (printIter) {
      cat("Iteration: ", format(iiter, width = nchar(iterlim)), "      Objective value: ", out$value, "\n")
    }
    
    if (accept) {
      
      if (minimize)
        g_boundary <- c(upper[which(out$gradient[upper] < 0)], lower[which(out$gradient[lower] > 0)])
      if (!minimize)
        g_boundary <- c(upper[which(out$gradient[upper] > 0)], lower[which(out$gradient[lower] < 0)])
      
      g_noboundary <- setdiff(1:d, g_boundary)
      n_boundary <- length(g_boundary)
      if (n_boundary > 0) {
        out$gradient <- out$gradient[-g_boundary]
        out$hessian <- out$hessian[-g_boundary, , drop = FALSE]
        out$hessian <- out$hessian[, -g_boundary, drop = FALSE]
      }
      
      
      B <- out$hessian
      g <- out$gradient
      f <- out$value
      out.value.save <- f
      if (rescale) {
        B <- B/outer(parscale, parscale)
        g <- g/parscale
      }
      if (!minimize) {
        B <- (-B)
        g <- (-g)
        f <- (-f)
      }
      eout <- eigen(B, symmetric = TRUE)
      gq <- as.numeric(t(eout$vectors) %*% g)
    }
    is.newton <- FALSE
    if (all(eout$values > 0)) {
      ptry <- as.numeric(-eout$vectors %*% (gq/eout$values))
      if (norm(ptry) <= r) 
        is.newton <- TRUE
    }
    if (!is.newton) {
      lambda.min <- min(eout$values)
      beta <- eout$values - lambda.min
      imin <- beta == 0
      C1 <- sum((gq/beta)[!imin]^2)
      C2 <- sum(gq[imin]^2)
      C3 <- sum(gq^2)
      if (C2 > 0 || C1 > r^2) {
        is.easy <- TRUE
        is.hard <- (C2 == 0)
        beta.dn <- sqrt(C2)/r
        beta.up <- sqrt(C3)/r
        fred <- function(beep) {
          if (beep == 0) {
            if (C2 > 0) 
              return(-1/r)
            else return(sqrt(1/C1) - 1/r)
          }
          return(sqrt(1/sum((gq/(beta + beep))^2)) - 
                   1/r)
        }
        if (fred(beta.up) <= 0) {
          uout <- list(root = beta.up)
        }
        else if (fred(beta.dn) >= 0) {
          uout <- list(root = beta.dn)
        }
        else {
          uout <- stats::uniroot(fred, c(beta.dn, beta.up))
        }
        wtry <- gq/(beta + uout$root)
        ptry <- as.numeric(-eout$vectors %*% wtry)
      }
      else {
        is.hard <- TRUE
        is.easy <- FALSE
        wtry <- gq/beta
        wtry[imin] <- 0
        ptry <- as.numeric(-eout$vectors %*% wtry)
        utry <- sqrt(r^2 - sum(ptry^2))
        if (utry > 0) {
          vtry <- eout$vectors[, imin, drop = FALSE]
          vtry <- vtry[, 1]
          ptry <- ptry + utry * vtry
        }
      }
    }
    preddiff <- sum(ptry * (g + as.numeric(B %*% ptry)/2))
    if (rescale) {
      theta.try <- theta
      theta.try[g_noboundary] <- theta.try[g_noboundary] + ptry/parscale[g_noboundary]
    }
    else {
      theta.try <- theta
      theta.try[g_noboundary] <- theta.try[g_noboundary] + ptry
    }
    upper <- which(!(theta.try < parupper))
    lower <- which(!(theta.try > parlower))
    theta.try[upper] <- parupper[upper]
    theta.try[lower] <- parlower[lower]
    
    out <- try(objfun(theta.try, ...))
    if (inherits(out, "try-error")) 
      break
    checks <- try(check.objfun.output(out, minimize, d))
    if (inherits(checks, "try-error")){
      out <- c(as.list(out), error = checks)
      break
    }
    
    ftry <- out$value
    if (!minimize) 
      ftry <- (-ftry)
    rho <- (ftry - f)/preddiff
    if (ftry < Inf) {
      is.terminate <- abs(ftry - f) < fterm || abs(preddiff) < 
        mterm
    }
    else {
      is.terminate <- FALSE
      rho <- (-Inf)
    }
    if (is.terminate) {
      if (ftry < f) {
        accept <- TRUE
        theta <- theta.try
      }
    }
    else {
      if (rho < 1/4) {
        accept <- FALSE
        r <- r/4
      }
      else {
        accept <- TRUE
        theta <- theta.try
        if (rho > 3/4 && (!is.newton)) 
          r <- min(2 * r, rmax)
      }
    }
    
    if (blather) {
      theta.try.blather <- rbind(theta.try.blather, theta.try)
      val.try.blather <- c(val.try.blather, out$value)
      accept.blather <- c(accept.blather, accept)
      preddiff.blather <- c(preddiff.blather, preddiff)
      stepnorm.blather <- c(stepnorm.blather, norm(ptry))
      if (is.newton) {
        mytype <- "Newton"
      }
      else {
        if (is.hard) {
          if (is.easy) {
            mytype <- "hard-easy"
          }
          else {
            mytype <- "hard-hard"
          }
        }
        else {
          mytype <- "easy-easy"
        }
      }
      type.blather <- c(type.blather, mytype)
      rho.blather <- c(rho.blather, rho)
    }
    if (is.terminate) 
      break
  }
  
  # Finalize ----
  if (inherits(out, "try-error")) {
    # error of last iteration
    out <- c(error = out, list(argument = theta.try), converged = FALSE)
  } else if (!is.null(out$error)){
    # error of last iteration
    out <- c(out, list(argument = theta.try), converged = FALSE)
  } else {
    out <- try(objfun(theta, ...))
    if (inherits(out, "try-error")) {
      # error of final argument
      out <- list(error = out, argument = theta)
      warning("error in last call to objfun")
    }
    else {
      # error of final argument
      checks <- try(check.objfun.output(out, minimize, d))
      if (inherits(checks, "try-error"))
        out <- c(as.list(out), error = checks, list(argument = theta))
    }
    out$argument <- theta
    out$converged <- is.terminate
  }
  out$iterations <- iiter
  if (blather) {
    dimnames(theta.blather) <- NULL
    out$argpath <- theta.blather
    dimnames(theta.try.blather) <- NULL
    out$argtry <- theta.try.blather
    out$steptype <- type.blather
    out$accept <- accept.blather
    out$r <- r.blather
    out$rho <- rho.blather
    out$valpath <- val.blather
    out$valtry <- val.try.blather
    if (!minimize) 
      preddiff.blather <- (-preddiff.blather)
    out$preddiff <- preddiff.blather
    out$stepnorm <- stepnorm.blather
  }
  return(out)
}


check.objfun.output <- function(obj, minimize, dimen)
{
  if (! is.list(obj))
    stop("objfun returned object that is not a list")
  foo <- obj$value
  if (is.null(foo))
    stop("objfun returned list that does not have a component 'value'")
  if (! is.numeric(foo))
    stop("objfun returned value that is not numeric")
  if (length(foo) != 1)
    stop("objfun returned value that is not scalar")
  if (is.na(foo) || is.nan(foo))
    stop("objfun returned value that is NA or NaN")
  if (minimize && foo == (-Inf))
    stop("objfun returned -Inf value in minimization")
  if ((! minimize) && foo == Inf)
    stop("objfun returned +Inf value in maximization")
  if (is.finite(foo)) {
    bar <- obj$gradient
    if (is.null(bar))
      stop("objfun returned list without component 'gradient' when value is finite")
    if (! is.numeric(bar))
      stop("objfun returned gradient that is not numeric")
    if (length(bar) != dimen)
      stop(paste("objfun returned gradient that is not vector of length", dimen))
    if (! all(is.finite(bar)))
      stop("objfun returned gradient not having all elements finite")
    baz <- obj$hessian
    if (is.null(baz))
      stop("objfun returned list without component 'hessian' when value is finite")
    if (! is.numeric(baz))
      stop("objfun returned hessian that is not numeric")
    if (! is.matrix(baz))
      stop("objfun returned hessian that is not matrix")
    if (! all(dim(baz) == dimen))
      stop(paste("objfun returned hessian that is not", dimen, "by", dimen, "matrix"))
    if (! all(is.finite(baz)))
      stop("objfun returned hessian not having all elements finite")
  }
  return(TRUE)
}


