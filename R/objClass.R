## Methods for class "objfn" -----------------------------------------------



## Class "objlist" and its constructors ------------------------------------



#' Generate objective list from numeric vector
#' 
#' @param p Named numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.objlist(p)
#' @export
as.objlist <- function(p) {
  
  objlist(value = 0,
      gradient = structure(rep(0, length(p)), names = names(p)),
      hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
}


#' L2 norm between data and model prediction
#' 
#' @description For parameter estimation and optimization, an objective function
#' is needed. \code{normL2} returns an objective function for the L2 norm of
#' data and model prediction. The resulting objective function can be used for
#' optimization with the trust optimizer, see \link{mstrust}.
#' @param data object of class \link{datalist}
#' @param x object of class \link{prdfn}
#' @param errmodel object of class \link{obsfn}
#' @param times numeric vector, the time points where the prediction function is to be
#' evaluated. If NULL, time points are extacted from the datalist. If the prediction
#' function makes use of events, \code{times} should be set by hand.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @return Object of class \code{obsfn}, i.e. a function 
#' \code{obj(..., fixed, deriv, conditions, env)} that returns an objective list,
#' \link{objlist}.
#' @details Objective functions can be combined by the "+" operator, see \link{sumobjfn}.
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL, attr.name = "data") {

  if (is.null(times)) timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time))))) else timesD <- times

  data.conditions <- names(data)
  
  controls <- list(times = timesD, attr.name = attr.name)
  
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = data.conditions, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    template <- objlist(
      value = 0,
      gradient = structure(rep(0, length(pouter)), names = names(pouter)),
      hessian = matrix(0, nrow = length(pouter), ncol = length(pouter), dimnames = list(names(pouter), names(pouter)))
    )
   
    # Import from controls
    timesD <- controls$times
    attr.name <- controls$attr.name
    
    # Create new environment if necessary
    if (is.null(env)) env <- new.env()
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) {
      err <- NULL
      if (!is.null(errmodel)) {
        err <- errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
        mywrss <- nll(res(data[[cn]], prediction[[cn]], err[[cn]]))
      } else {
        mywrss <- wrss(res(data[[cn]], prediction[[cn]]))  
      }
      available <- intersect(names(pouter), names(mywrss$gradient))
      result <- template
      result$value <- mywrss$value
      result$gradient[available] <- mywrss$gradient[available]
      result$hessian[available, available] <- mywrss$hessian[available, available]
      return(result)
    })
    out.data <- Reduce("+", out.data)
    
    # Combine contributions and attach attributes
    out <- out.data
    attr(out, controls$attr.name) <- out.data$value
    
    if (!is.null(env)) {
      assign("prediction", prediction, envir = env)
    }
    
    attr(out, "env") <- env
    return(out)
      

  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(x, "parameters")
  return(myfn)

}


#' Soft L2 constraint on parameters
#' 
#' @param mu named numeric, the prior values
#' @param sigma named numeric of length of mu or numeric of length one.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return object of class \code{objfn}
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{p-\mu}{\sigma}\right)^2}{(p-mu)^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' mu <- c(A = 0, B = 0)
#' sigma <- c(A = 0.1, B = 1)
#' myfn <- constraintL2(mu, sigma)
#' myfn(pars = c(A = 1, B = -1))
#' @export
constraintL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {

  ## Augment sigma if length = 1
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  
  
  controls <- list(mu = mu, sigma = sigma, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Import from controls
    mu <- controls$mu
    sigma <- controls$sigma
    attr.name <- controls$attr.name
    
    # pouter can be a list (if result from a parameter transformation)
    # In this case match with conditions and evaluate only those
    # If there is no overlap, return NULL
    # If pouter is not a list, evaluate the constraint function 
    # for this pouter.
    
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      defined <- ifelse(is.null(condition), TRUE, condition %in% conditions)
      
      if (length(available) == 0 | !defined) return()
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
      par.fixed <- intersect(names(mu), names(fixed))
      sumOfFixed <- 0
      if (!is.null(par.fixed)) sumOfFixed <- sum(0.5*((fixed[par.fixed] - mu[par.fixed])/sigma[par.fixed]) ^ 2)
      
      # Compute prior value and derivatives
      par <- intersect(names(mu), names(p))
      
      val <- sum((0.5*((p[par] - mu[par])/sigma[par]) ^ 2)) + sumOfFixed
      gr <- rep(0, length(p)); names(gr) <- names(p)
      gr[par] <- ((p[par] - mu[par])/(sigma[par] ^ 2))
      
      hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
      diag(hs)[par] <- 1/sigma[par] ^ 2
      
      dP <- attr(p, "deriv")
      if (!is.null(dP)) {
        gr <- as.vector(gr %*% dP); names(gr) <- colnames(dP)
        hs <- t(dP) %*% hs %*% dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
      }
      
      objlist(value = 2*val, gradient = 2*gr, hessian = 2*hs)
      
      
    })
    
    out <- Reduce("+", outlist)
    attr(out, controls$attr.name) <- out$value
    attr(out, "env") <- env
    return(out)
    
    
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  return(myfn)
 
  
}




#' L2 objective function for validation data point
#' 
#' @param name character, the name of the prediction, e.g. a state name.
#' @param time numeric, the time-point associated to the prediction
#' @param value character, the name of the parameter which contains the
#' prediction value.
#' @param sigma numeric, the uncertainty of the introduced test data point
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the prediction is made.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintL2}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' prediction <- list(a = matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A"))))
#' derivs <- matrix(c(0, 1, 0.1), nrow = 1, dimnames = list(NULL, c("time", "A.A", "A.k1")))
#' attr(prediction$a, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' 
#' vali <- datapointL2(name = "A", time = 0, value = "newpoint", sigma = 1, condition = "a")
#' vali(pars = c(p0, newpoint = 1), env = .GlobalEnv)
#' @export
datapointL2 <- function(name, time, value, sigma = 1, attr.name = "validation", condition) {
  
  
  controls <- list(
    mu = structure(name, names = value)[1], # Only one data point is allowed
    time = time[1],
    sigma = sigma[1],
    attr.name = attr.name
  )
  
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = NULL, env = NULL) {
  
    # Import from controls
    mu <- controls$mu
    time <- controls$time
    sigma <- controls$sigma
    attr.name <- controls$attr.name
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    if (!is.null(env)) {
      prediction <- as.list(env)$prediction 
    } else {
      stop("No prediction available. Use the argument env to pass an environment that contains the prediction.")
    }
  
    # Return result only when the data point condition overlaps with the evaluated conditions
    if (!is.null(conditions) && !condition %in% conditions) return()
    
    # Divide parameter into data point and rest
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    
    
    # Get predictions and derivatives at time point
    time.index <- which(prediction[[condition]][,"time"] == time)
    withDeriv <- !is.null(attr(prediction[[condition]], "deriv"))
    pred <- prediction[[condition]][time.index, ]
    deriv <- NULL
    if (withDeriv)
      deriv <- attr(prediction[[condition]], "deriv")[time.index, ]
    
    # Reduce to name = mu
    pred <- pred[mu]
    if (withDeriv) {
      mu.para <- intersect(paste(mu, parapar, sep = "."), names(deriv))
      deriv <- deriv[mu.para]
    }
    
    # Compute prior value and derivatives
    res <- as.numeric(pred - c(fixed, pouter)[names(mu)])
    val <- as.numeric((res/sigma) ^ 2)
    gr <- NULL
    hs <- NULL
    
    if (withDeriv) {
      dres.dp <- structure(rep(0, length(pouter)), names = names(pouter))
      if (length(parapar) > 0) dres.dp[parapar] <- as.numeric(deriv)
      if (length(datapar) > 0) dres.dp[datapar] <- -1
      gr <- 2*res*dres.dp/sigma ^ 2
      hs <- 2*outer(dres.dp, dres.dp, "*")/sigma ^ 2; colnames(hs) <- rownames(hs) <- names(pouter)
    }
    
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, controls$attr.name) <- out$value
    
    attr(out, "env") <- env
    
    return(out)
    
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  
  
  return(myfn)
  
      
}

#' L2 objective function for prior value
#' 
#' @description As a prior function, it returns derivatives with respect to
#' the penalty parameter in addition to parameter derivatives.
#' 
#' @param mu Named numeric, the prior values
#' @param lambda Character of length one. The name of the penalty paramter in \code{p}.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{e^{\lambda} \| p-\mu \|^2}{exp(lambda)*sum((p-mu)^2)}
#' and its derivatives with respect to p and lambda.
#' @examples
#' p <- c(A = 1, B = 2, C = 3, lambda = 0)
#' mu <- c(A = 0, B = 0)
#' obj <- priorL2(mu = mu, lambda = "lambda")
#' obj(pars = p + rnorm(length(p), 0, .1))
#' @export
priorL2 <- function(mu, lambda = "lambda", attr.name = "prior", condition = NULL) {
  
  
  controls <- list(mu = mu, lambda = lambda, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
   
    # Import from controls 
    mu <- controls$mu
    lambda <- controls$lambda
    attr.name <- controls$attr.name
    
    # pouter can be a list (if result from a parameter transformation)
    # In this case match with conditions and evaluate only those
    # If there is no overlap, return NULL
    # If pouter is not a list, evaluate the constraint function 
    # for this pouter.
    
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      defined <- ifelse(is.null(condition), TRUE, condition %in% conditions)
      
      if (length(available) == 0 | !defined) return()
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      
      ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
      par.fixed <- intersect(names(mu), names(fixed))
      sumOfFixed <- 0
      if (!is.null(par.fixed)) sumOfFixed <- sum(exp(c(fixed, p)[lambda])*(fixed[par.fixed] - mu[par.fixed]) ^ 2)
      
      # Compute prior value and derivatives
      par <- intersect(names(mu), names(p))
      par0 <- setdiff(par, lambda)
      
      val <- sum(exp(c(fixed, p)[lambda]) * (p[par] - mu[par]) ^ 2) + sumOfFixed
      gr <- rep(0, length(p)); names(gr) <- names(p)
      gr[par] <- 2*exp(c(fixed, p)[lambda])*(p[par] - mu[par])
      if (lambda %in% names(p)) {
        gr[lambda] <- sum(exp(c(fixed, p)[lambda]) * (p[par0] - mu[par0]) ^ 2) + 
          sum(exp(c(fixed, p)[lambda]) * (fixed[par.fixed] - mu[par.fixed]) ^ 2)
      }
      
      hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
      diag(hs)[par] <- 2*exp(c(fixed, p)[lambda])
      if (lambda %in% names(p)) {
        hs[lambda, lambda] <- gr[lambda] 
        hs[lambda, par0] <- hs[par0, lambda] <- gr[par0]
      }
      
      dP <- attr(p, "deriv")
      if (!is.null(dP)) {
        gr <- as.vector(gr %*% dP); names(gr) <- colnames(dP)
        hs <- t(dP) %*% hs %*% dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
      }
      
      objlist(value = val, gradient = gr, hessian = hs)
      
    })
    
    out <- Reduce("+", outlist)
    attr(out, controls$attr.name) <- out$value
    attr(out, "env") <- env
    
    return(out)
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  return(myfn)
  
}


#' Compute the weighted residual sum of squares
#' 
#' @param nout data.frame (result of \link{res}) or object of class \link{objframe}.
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
#' @export
wrss <- function(nout) {
  
  
  
  obj <- sum(nout$weighted.residual^2)
  grad <- NULL
  hessian <- NULL
  
  if (!is.null(attr(nout, "deriv"))) {
    
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
  
    sens <- as.matrix(attr(nout, "deriv")[, -(1:2), drop = FALSE])
    grad <- as.vector(2*matrix(nout$residual/nout$sigma^2, nrow = 1) %*% sens)
    names(grad) <- colnames(sens)
    hessian <- 2*t(sens/nout$sigma) %*% (sens/nout$sigma)
    
  }
  
  
  objlist(value = obj, gradient = grad, hessian = hessian)
  
}
#' Compute the negative log-likelihood
#' 
#' @param nout data.frame (result of \link{res}) or object of class \link{objframe}.
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
#' @export
nll <- function(nout) {
  
  obj <- sum(nout$weighted.residual^2) + 2 * sum(log(nout$sigma))
  grad <- NULL
  hessian <- NULL
  
  
  if (!is.null(attr(nout, "deriv")) & !is.null(attr(nout, "deriv.err"))) {
    
    
    sens <- as.matrix(attr(nout, "deriv")[, -(1:2), drop = FALSE])
    sens.err <- as.matrix(attr(nout, "deriv.err")[, -(1:2), drop = FALSE])
    
    grad <- as.vector(2*matrix(nout$residual/nout$sigma^2, nrow = 1) %*% sens -
                        2*matrix(nout$residual^2/nout$sigma^3, nrow = 1) %*% sens.err +
                        2*matrix(1/nout$sigma, nrow = 1) %*% sens.err)
    names(grad) <- colnames(sens)
    
    hessian <- 2 * t(sens/nout$sigma - sens.err*nout$residual/nout$sigma^2) %*% (sens/nout$sigma - sens.err*nout$residual/nout$sigma^2) # - 2 * t(sens.err/nout$sigma) %*% (sens.err/nout$sigma)
    
  }
  
  
  objlist(value = obj, gradient = grad, hessian = hessian)
  
}


## Methods for class objlist ------------------------------------------------

#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
#' @aliases sumobjlist
#' @export "+.objlist"
#' @export
"+.objlist" <- function(out1, out2) {
  
  if (is.null(out1)) return(out2)
  if (is.null(out2)) return(out1)

  
  
  allnames <- c(names(out1), names(out2))
  what <- allnames[duplicated(allnames)]
  what.names <- what
  if (is.null(what)) {
    what <- 1:min(c(length(out1), length(out2)))
    what.names <- NULL
  }
  
  out12 <- lapply(what, function(w) {
    sub1 <- out1[[w]]
    sub2 <- out2[[w]]
    n <- names(sub1)
    dn <- dimnames(sub1)
    if (!is.null(n) && !is.null(sub1) %% !is.null(sub2)) {
      #print("case1: sum of vectors")
      sub1[n] + sub2[n]
    } else if (!is.null(dn) && !is.null(sub1) && !is.null(sub2)) {
      #print("case2: sum of matrices")
      matrix(sub1[dn[[1]], dn[[2]]] + sub2[dn[[1]], dn[[2]]], 
             length(dn[[1]]), length(dn[[2]]), dimnames = list(dn[[1]], dn[[2]]))
    } else if (!is.null(sub1) && !is.null(sub2)) {
      #print("case3: sum of scalars")
      sub1 + sub2
    } else {
      #print("case4")
      NULL
    }
  })
  names(out12) <- what.names
  
  # Summation of numeric attributes 
  out1.attributes <- attributes(out1)[sapply(attributes(out1), is.numeric)]
  out2.attributes <- attributes(out2)[sapply(attributes(out2), is.numeric)]
  attr.names <- union(names(out1.attributes), names(out2.attributes))
  out12.attributes <- lapply(attr.names, function(n) {
    x1 <- ifelse(is.null(out1.attributes[[n]]), 0, out1.attributes[[n]])
    x2 <- ifelse(is.null(out2.attributes[[n]]), 0, out2.attributes[[n]])
    x1 + x2
  })
  attributes(out12)[attr.names] <- out12.attributes
  
  
  
  class(out12) <- "objlist"
  
  return(out12)
}


#' @export
print.objfn <- function(x, ...) {
 
  parameters <- attr(x, "parameters")
  
  cat("Objective function:\n")
  str(args(x))
  cat("\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  
}

