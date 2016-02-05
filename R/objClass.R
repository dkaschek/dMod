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
#' as.obj(p)
#' @export
as.objlist <- function(p) {
  
  objlist(value = 0,
      gradient = structure(rep(0, length(p)), names = names(p)),
      hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
}


#' Compute a differentiable box prior
#' 
#' @param p Named numeric, the parameter value
#' @param mu Named numeric, the prior values, means of boxes
#' @param sigma Named numeric, half box width
#' @param k Named numeric, shape of box; if 0 a quadratic prior is obtained, the higher k the more box shape, gradient at border of the box (-sigma, sigma) is equal to sigma*k
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value but not to gradient and Hessian)
#' @return list with entries: value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric). Object of class \code{objlist}.
#' @export
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, fixed=NULL) {
  
  kmin <- 1e-5
  
  ## Augment sigma if length = 1
  if(length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu)) 
  ## Augment k if length = 1
  if(length(k) == 1) 
    k <- structure(rep(k, length(mu)), names = names(mu))
  
  k <- sapply(k, function(ki){
    if(ki < kmin){
      kmin
    } else ki
  })
  
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(0.5*(exp(k[par.fixed]*((fixed[par.fixed] - mu[par.fixed])/sigma[par.fixed])^2)-1)/(exp(k[par.fixed])-1))
  
  
  par <- intersect(names(mu), names(p))
  t <- p[par]
  mu <- mu[par]
  s <- sigma[par]
  k <- k[par]
  
  # Compute prior value and derivatives 
  
  gr <- rep(0, length(t)); names(gr) <- names(t)
  hs <- matrix(0, length(t), length(t), dimnames = list(names(t), names(t)))
  
  val <- sum(0.5*(exp(k*((t-mu)/s)^2)-1)/(exp(k)-1)) + sumOfFixed
  gr <- (k*(t-mu)/(s^2)*exp(k*((t-mu)/s)^2)/(exp(k)-1))
  diag(hs)[par] <- k/(s*s)*exp(k*((t-mu)/s)^2)/(exp(k)-1)*(1+2*k*(t-mu)/(s^2))
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  objlist(value=val,gradient=gr,hessian=hs)
  
}



#' @export
normL2 <- function(data, x, times = NULL) {

  if (is.null(times)) timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time))))) else timesD <- times

  myfn <- function(pouter, fixed = NULL, deriv=TRUE, conditions = names(data), env = NULL) {
    
    # Evaluate the code in the local environment or in env
    # If evaluated in env, arguments have to be transferred
    # to this environment.
    if (is.null(env)) {
      env <- environment()
    } else {
      assign("pouter", pouter, envir = env)
      assign("fixed", fixed, envir = env)
      assign("deriv", deriv, envir = env)
      assign("conditions", conditions, envir = env)
      assign("timesD", timesD, envir = env)
      assign("x", x, envir = env)
      assign("data", data, envir = env) 
      
    }
    
    myexpression <- expression({
    
      prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)
      
      # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
      out.data <- lapply(conditions, function(cn) wrss(res(data[[cn]], prediction[[cn]])))
      out.data <- Reduce("+", out.data)
      
      # Combine contributions and attach attributes
      out <- out.data
      attr(out, "data") <- out.data$value
      attr(out, "env") <- env
      
      return(out)
      
    })
    
    eval(myexpression, env, parent.frame())
    

  }
  class(myfn) <- c("objfn", "fn")
  return(myfn)

}


#' Soft L2 constraint on parameters
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param sigma Named numeric of length of mu or numeric of length one.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{\frac{1}{2}\left(\frac{p-\mu}{\sigma}\right)^2}{0.5*(p-mu)^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' p <- c(A = 1, B = 2, C = 3)
#' mu <- c(A = 0, B = 0)
#' sigma <- c(A = 0.1, B = 1)
#' constraintL2(p, mu, sigma)
#' @export
constraintL2 <- function(mu, sigma = 1, condition = NULL) {

  ## Augment sigma if length = 1
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  
  
  myfn <- function(pouter, fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    # Evaluate the code in the local environment or in env
    # If evaluated in env, arguments have to be transferred
    # to this environment.
    if (is.null(env)) {
      env <- environment()
    } else {
      assign("pouter", pouter, envir = env)
      assign("fixed", fixed, envir = env)
      assign("deriv", deriv, envir = env)
      assign("condition", condition, envir = env)
      assign("conditions", conditions, envir = env)
      assign("mu", mu, envir = env)
      assign("sigma", sigma, envir = env)
    }
    
    myexpression <- expression({
      
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
        
        objlist(value = val, gradient = gr, hessian = hs)
        
        
      })
      
      out <- Reduce("+", outlist)
      attr(out, "prior") <- out$value
      attr(out, "env") <- env
      return(out)
      
    })
  
    eval(myexpression, env, parent.frame())
    
    
  }
  class(myfn) <- c("objfn", "fn")
  return(myfn)
 
  
}




#' L2 objective function for validation data point
#' 
#' @param p Namec numeric, the parameter values
#' @param prediction Matrix with first column "time" and one column per predicted state. Can have
#' an attribute \code{deriv}, the matrix of sensitivities. If present, derivatives of the objective
#' function with respect to the parameters are returned.
#' @param mu Named character of length one. Has the structure \code{mu = c(parname = statename)}, where
#' \code{statename} is one of the column names of \code{prediction} and \code{parname} is one of the
#' names of \code{p}, allowing to treat the validation data point as a parameter.
#' @param time Numeric of length one. An existing time point in \code{prediction}.
#' @param sigma Numeric of length one. The uncertainty assumed for the validation data point.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintL2}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' \dontrun{
#' prediction <- matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A")))
#' derivs <- matrix(c(0, 1, 0.1), nrow = 1, dimnames = list(NULL, c("time", "A.A", "A.k1")))
#' attr(prediction, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' mu <- c(newpoint = "A")
#' timepoint <- 0
#' 
#' datapointL2(p = c(p, newpoint = 2), prediction, mu, timepoint)
#' datapointL2(p = c(p, newpoint = 1), prediction, mu, timepoint)
#' datapointL2(p = c(p, newpoint = 0), prediction, mu, timepoint)
#' }
#' @export
datapointL2 <- function(name, time, value, sigma = 1, condition) {
  
  # Only one data point is allowed
  mu <- structure(name, names = value)
  mu <- mu[1]; time <- time[1]; sigma <- sigma[1]
  
  
  myfn <- function(pouter, fixed = NULL, deriv=TRUE, conditions = NULL, env = NULL) {
  
    # Evaluate the code in the local environment or in env
    # If evaluated in env, arguments have to be transferred
    # to this environment.
    if (is.null(env)) {
      env <- environment()
    } else {
      assign("pouter", pouter, envir = env)
      assign("fixed", fixed, envir = env)
      assign("deriv", deriv, envir = env)
      assign("conditions", conditions, envir = env)
      assign("mu", mu, envir = env)
      assign("time", time, envir = env)
      assign("sigma", sigma, envir = env)
      assign("condition", condition, envir = env)
    }

    myexpression <- expression({    
  
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
      attr(out, "validation") <- out$value
      
      attr(out, "env") <- env
      
      return(out)
      
    })
    
    eval(myexpression, env, parent.frame())
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  
  return(myfn)
  
      
}

#' L2 objective function for prior value
#' 
#' @description As a prior function, it returns derivatives with respect to
#' the penalty parameter in addition to parameter derivatives.
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param lambda Character of length one. The name of the penalty paramter in \code{p}.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintExp2}
#' @details Computes the constraint value 
#' \deqn{e^{\lambda} \| p-\mu \|^2}{exp(lambda)*sum((p-mu)^2)}
#' and its derivatives with respect to p and lambda.
#' @examples
#' p <- c(A = 1, B = 2, C = 3, lambda = 0)
#' mu <- c(A = 0, B = 0)
#' priorL2(p, mu, lambda = "lambda")
#' @export
priorL2 <- function(mu, lambda = "lambda", condition = NULL) {
  
  
  myfn <- function(pouter, fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    # Evaluate the code in the local environment or in env
    # If evaluated in env, arguments have to be transferred
    # to this environment.
    if (is.null(env)) {
      env <- environment()
    } else {
      assign("pouter", pouter, envir = env)
      assign("fixed", fixed, envir = env)
      assign("deriv", deriv, envir = env)
      assign("condition", condition, envir = env)
      assign("conditions", conditions, envir = env)
      assign("mu", mu, envir = env)
      assign("lambda", lambda, envir = env)
    }
    
    myexpression <- expression({
      
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
      attr(out, "prior") <- out$value
      attr(out, "env") <- env
      
      return(out)
      
    })
    
    eval(myexpression, env, parent.frame())
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
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
  
  if(!is.null(attr(nout, "deriv"))) {
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
  
    sens <- as.matrix(attr(nout, "deriv")[,-(1:2)])
    grad <- as.vector(2*matrix(nout$residual/nout$sigma^2, nrow=1)%*%sens)
    names(grad) <- colnames(sens)
    hessian <- 2*t(sens/nout$sigma)%*%(sens/nout$sigma)
    
    
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
#' @aliases summation
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

