
norm <- function(x) sqrt(sum(x^2))

#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
"+.obj" <- function(out1, out2) {
  
  allnames <- c(names(out1), names(out2))
  what <- allnames[duplicated(allnames)]
  what.names <- what
  if(is.null(what)) {
    what <- 1:min(c(length(out1), length(out2)))
    what.names <- NULL
  }
  
  out12 <- lapply(what, function(w) {
    sub1 <- out1[[w]]
    sub2 <- out2[[w]]
    n <- names(sub1)
    dn <- dimnames(sub1)
    if(!is.null(n) && !is.null(sub1) %% !is.null(sub2)) {
      #print("case1")
      sub1[n] + sub2[n]
    } else if(!is.null(dn) && !is.null(sub1) && !is.null(sub2)) {
      #print("case2")
      sub1[dn[[1]], dn[[2]]] + sub2[dn[[1]], dn[[2]]]
    } else if(!is.null(sub1) && !is.null(sub2)) {
      #print("case3")
      sub1 + sub2
    } else {
      #print("case4")
      NULL
    }
  })
  names(out12) <- what.names
  
  class(out12) <- c("obj", "list")

  return(out12)
}


########## REFERENCES ##########
#####
##### Fletcher, R. (1987)
##### Practical Methods of Optimization, second edition.
##### John Wiley, Chichester.
#####
##### Nocedal, J. and Wright, S. J. (1999)
##### Numerical Optimization.
##### Springer-Verlag, New York.
#####
##### See Section 5.1 of Fletcher
##### See Section 4.2 of Nocedal and Wright
#####
################################

########## COMMENT ##########
##### Our method using one eigendecomposition per iteration is not fastest.
##### Both books recommend using multiple Cholesky decompositions instead.
##### But the eigendecomposition method is simpler to program, easier to
##### understand (which is why both books use it for their theoretical
##### explanation), and hopefully more bulletproof.
#####
##### Our idea for this comes from the way mvrnorm in the MASS package also
##### uses eigendecomposition rather than Cholesky -- also because bulletproof
##### is better than fast.
#############################

trustL1 <- function(objfun, parinit, mu = 0*parinit, lambda = 1, rinit, rmax, parscale,
    iterlim = 100, fterm = sqrt(.Machine$double.eps),
    mterm = sqrt(.Machine$double.eps),
    minimize = TRUE, blather = FALSE, ...)
{
    if (! is.numeric(parinit))
       stop("parinit not numeric")
    if (! all(is.finite(parinit)))
       stop("parinit not all finite")
    d <- length(parinit)
    if (missing(parscale)) {
        rescale <- FALSE
    } else {
        rescale <- TRUE
        if (length(parscale) != d)
           stop("parscale and parinit not same length")
        if (! all(parscale > 0))
           stop("parscale not all positive")
        if (! all(is.finite(parscale) & is.finite(1 / parscale)))
           stop("parscale or 1 / parscale not all finite")
    }
    if (! is.logical(minimize))
       stop("minimize not logical")

    r <- rinit
    theta <- parinit
    out <- try(objfun(theta, ...))
    grad0 <- out$gradient
    outL1 <- try(constraintL1(theta, names(mu), mu, lambda))
    gradL1 <- outL1$gradient
    
    if (inherits(out, "try-error")) {
        warning("error in first call to objfun")
        return(list(error = out, argument = theta, converged = FALSE,
            iterations = 0))
    }
    
    out <- out + outL1
    check.objfun.output(out, minimize, d)
    if (! is.finite(out$value))
        stop("parinit not feasible")
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

    for (iiter in 1:iterlim) {

        if (blather) {
            theta.blather <- rbind(theta.blather, theta)
            r.blather <- c(r.blather, r)
            if (accept)
                val.blather <- c(val.blather, out$value)
            else
                val.blather <- c(val.blather, out.value.save)
        }

        if (accept) {
            B <- out$hessian
            g <- out$gradient
            f <- out$value
            out.value.save <- f
            if (rescale) { 
                B <- B / outer(parscale, parscale)
                g <- g / parscale
            }
            if (! minimize) {
                B <- (- B)
                g <- (- g)
                f <- (- f)
            }
            eout <- eigen(B, symmetric = TRUE)
            gq <- as.numeric(t(eout$vectors) %*% g)
        }

        ########## solve trust region subproblem ##########

        ##### try for Newton #####
        is.newton <- FALSE
        if (all(eout$values > 0)) {
            ptry <- as.numeric(- eout$vectors %*% (gq / eout$values))
            if (norm(ptry) <= r)
                is.newton <- TRUE
        }

        ##### non-Newton #####
        if (! is.newton) {
            lambda.min <- min(eout$values)
            beta <- eout$values - lambda.min
            imin <- beta == 0
            C1 <- sum((gq / beta)[! imin]^2)
            C2 <- sum(gq[imin]^2)
            C3 <- sum(gq^2)
            if (C2 > 0 || C1 > r^2) {
                is.easy <- TRUE
                is.hard <- (C2 == 0)
                ##### easy cases #####
                beta.dn <- sqrt(C2) / r
                beta.up <- sqrt(C3) / r
                fred <- function(beep) {
                    if (beep == 0) {
                        if (C2 > 0)
                            return(- 1 / r)
                        else
                            return(sqrt(1 / C1) - 1 / r)
                    }
                    return(sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r)
                }
                if (fred(beta.up) <= 0) {
                    uout <- list(root = beta.up)
                } else if (fred(beta.dn) >= 0) {
                    uout <- list(root = beta.dn)
                } else {
                    uout <- uniroot(fred, c(beta.dn, beta.up))
                }
                wtry <- gq / (beta + uout$root)
                ptry <- as.numeric(- eout$vectors %*% wtry)
            } else {
                is.hard <- TRUE
                is.easy <- FALSE
                ##### hard-hard case #####
                wtry <- gq / beta
                wtry[imin] <- 0
                ptry <- as.numeric(- eout$vectors %*% wtry)
                utry <- sqrt(r^2 - sum(ptry^2))
                if (utry > 0) {
                    vtry <- eout$vectors[ , imin, drop = FALSE]
                    vtry <- vtry[ , 1]
                    ptry <- ptry + utry * vtry
                }
            }
        }

        
        ########## predicted versus actual change ##########

        ## Check if step away from prior would lead to step-back
        names(ptry) <- names(theta) <- names(parinit)
        
        is.mu <- which(theta[names(mu)] == mu)
        subgradL1 <- lambda*sign(ptry[names(mu)][is.mu])
        subgrad <- grad0[names(mu)][is.mu]
        not.doStep <- sign(subgradL1 + subgrad)*sign(subgradL1) == 1
        ptry[names(mu)][not.doStep] <- 0
        
        ## Compute theta.try
        preddiff <- sum(ptry * (g + as.numeric(B %*% ptry) / 2))
        if (rescale) {
            theta.try <- theta + ptry / parscale
        } else {
            theta.try <- theta + ptry
        }
        
        ## Set on prior value if step-over
        chgsgn <- (theta[names(mu)]-mu)*(theta.try[names(mu)]-mu)
        theta.try[names(mu)][chgsgn < 0] <- mu[chgsgn < 0]
        
        
        out <- try(objfun(theta.try, ...))
        outL1 <- try(constraintL1(theta.try, mu, lambda))
        if (inherits(out, "try-error"))
            break
        out <- out + outL1
        check.objfun.output(out, minimize, d)
        ftry <- out$value
        if (! minimize)
            ftry <- (- ftry)
        rho <- (ftry - f) / preddiff

        ########## termination test ##########
        if (ftry < Inf) {
            is.terminate <- abs(ftry - f) < fterm || abs(preddiff) < mterm
        } else {
            is.terminate <- FALSE
            rho <- (- Inf)
        }

        ##### adjustments #####
        if (is.terminate) {
            if (ftry < f) {
                accept <- TRUE
                theta <- theta.try
            }
        } else {
            if (rho < 1 / 4) {
                accept <- FALSE
                r <- r / 4
            } else {
                accept <- TRUE
                theta <- theta.try
                if (rho > 3 / 4 && (! is.newton))
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
            } else {
                if (is.hard) {
                    if (is.easy) {
                        mytype <- "hard-easy"
                    } else {
                        mytype <- "hard-hard"
                    }
                } else {
                    mytype <- "easy-easy"
                }
            }
            type.blather <- c(type.blather, mytype)
            rho.blather <- c(rho.blather, rho)
        }

        if (is.terminate)
            break
    }

    if (inherits(out, "try-error")) {
        out <- list(error = out, argument = theta.try, converged = FALSE)
    } else {
        out <- try(objfun(theta, ...))
        outL1 <- try(constraintL1(theta.try, mu, lambda))
        if (inherits(out, "try-error")) {
          out <- list(error = out)
          warning("error in last call to objfun")
        } else {
          out <- out + outL1
          check.objfun.output(out, minimize, d)
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
        if (! minimize)
            preddiff.blather <- (- preddiff.blather)
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



constraintL1 <- function(p, mu, lambda = 1) {
  
  parameters <- names(mu)
  value <- sum(lambda*abs(p[parameters] - mu))
  
  gradient <- rep(0, length(p)); names(gradient) <- names(p)
  gradient[parameters][p[parameters] >  mu] <-  lambda
  gradient[parameters][p[parameters] <  mu] <- -lambda
  
  hessian <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
  diag(hessian)[parameters] <- lambda^2
  
  dP <- attr(p, "deriv") 
  if(!is.null(dP)) {
    gradient <- as.vector(gradient %*% dP)
    names(gradient) <- colnames(dP)
    hessian <- t(dP) %*% hessian %*% dP
    colnames(hessian) <- colnames(dP)
    rownames(hessian) <- colnames(dP)
  }
  
  out <- list(value = value, gradient = gradient, hessian = hessian)
  class(out) <- c("obj", "list")
  
  return(out)
  
  
}


#' Compute the L1 norm of residuals
#' 
#' @param nout data.frame (result of \link{res})
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
l1norm <- function(nout) {
  
  sign.res <- sign(nout$weighted.residual)
  obj <- sum(abs(nout$weighted.residual))
  grad <- NULL
  hessian <- NULL
  
  if(!is.null(attr(nout, "deriv"))) {
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
    
    sens <- as.matrix(attr(nout, "deriv")[,-(1:2)])
    grad <- apply(-sign.res*sens/nout$sigma, 2, sum)
    names(grad) <- colnames(sens)
    hessian <- t(sens/nout$sigma)%*%(sens/nout$sigma)
    
    
  }
  
  out <- list(value=obj, gradient=grad, hessian=hessian)
  class(out) <- c("obj", "list")
  
  return(out)
  
}

