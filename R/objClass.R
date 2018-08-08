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
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, fixed=NULL) {
  
  
  ##
  ## This function need to be extended according to constraintL2()
  ## The parameters sigma and k need to be replaced by more
  ## meaningful parameters.
  ##
  
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


#' L2 norm between data and model prediction
#' 
#' @description For parameter estimation and optimization, an objective function
#' is needed. \code{normL2} returns an objective function for the L2 norm of
#' data and model prediction. The resulting objective function can be used for
#' optimization with the trust optimizer, see \link{mstrust}.
#' @param data object of class \link{datalist}
#' @param x object of class \link{prdfn}
#' @param errmodel object of class \link{obsfn}. \code{errmodel} does not need to be defined for all conditions.
#' @param times numeric vector, additional time points where the prediction function is 
#' evaluated. If NULL, time points are extacted from the datalist solely. If the prediction
#' function makes use of events, hand over event \code{times} here.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param loq named numeric or single valued numeric. Limit of quantification.
#' @return Object of class \code{obsfn}, i.e. a function 
#' \code{obj(..., fixed, deriv, conditions, env)} that returns an objective list,
#' \link{objlist}.
#' @details Objective functions can be combined by the "+" operator, see \link{sumobjfn}.
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL, attr.name = "data", loq = -Inf) {

  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times)) timesD <- sort(union(times, timesD))

  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  if (!all(data.conditions %in% x.conditions)) 
    stop("The prediction function does not provide predictions for all conditions in the data.")
  e.conditions <- names(attr(errmodel, "mappings"))
  
  controls <- list(times = timesD, attr.name = attr.name, conditions = x.conditions, loq = loq)

  # might be necessary to "store" errmodel in the objective function (-> runbg)
  force(errmodel)  
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = controls$conditions, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Generate output template
    pars_out <- colnames(getDerivs(as.parvec(pouter)))
    template <- objlist(
      value = 0,
      gradient = structure(rep(0, length(pars_out)), names = pars_out),
      hessian = matrix(0, nrow = length(pars_out), ncol = length(pars_out), dimnames = list(pars_out, pars_out))
    )
   
    # Import from controls
    timesD <- controls$times
    attr.name <- controls$attr.name
    loq <- controls$loq
    
    # Create new environment if necessary
    if (is.null(env)) env <- new.env()
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) {
      err <- NULL
      if ((!is.null(errmodel) & is.null(e.conditions)) | (!is.null(e.conditions) && (cn %in% e.conditions))) {
        err <- errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
        mywrss <- nll(res(data[[cn]], prediction[[cn]], err[[cn]], loq))
      } else {
        mywrss <- wrss(res(data[[cn]], prediction[[cn]], NULL, loq))  
      }
      available <- intersect(pars_out, names(mywrss$gradient))
      result <- template
      result$value <- mywrss$value
      if (deriv) {
        result$gradient[available] <- mywrss$gradient[available]
        result$hessian[available, available] <- mywrss$hessian[available, available]  
      } else {
        result$gradient <- result$hessian <- NULL
      }
      
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
  attr(myfn, "modelname") <- modelname(x, errmodel)
  return(myfn)

}


#' Soft L2 constraint on parameters
#' 
#' @param mu named numeric, the prior values
#' @param sigma named numeric of length of mu or numeric of length one
#' or character of length of mu or character of length one
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return object of class \code{objfn}
#' @seealso \link{wrss}
#' @details If sigma is numeric, the function computes the constraint value 
#' \deqn{\left(\frac{p-\mu}{\sigma}\right)^2}{(p-mu)^2/sigma^2}
#' and its derivatives with respect to p. If sigma is a character, the 
#' function computes
#' \deqn{\left(\frac{p-\mu}{\sigma}\right)^2 + \log(\sigma^2)}{(p-mu)^2/sigma^2 + log(sigma^2)}
#' and its derivatives with respect to p and sigma. Sigma parameters being
#' passed to the function are ALWAYS assumed to be on a log scale, i.e. internally
#' sigma parameters are converted by \code{exp()}.
#' @examples
#' mu <- c(A = 0, B = 0)
#' sigma <- c(A = 0.1, B = 1)
#' myfn <- constraintL2(mu, sigma)
#' myfn(pars = c(A = 1, B = -1))
#' 
#' # Introduce sigma parameter but fix them (sigma parameters
#' # are assumed to be passed on log scale)
#' mu <- c(A = 0, B = 0)
#' sigma <- paste("sigma", names(mu), sep = "_")
#' myfn <- constraintL2(mu, sigma)
#' pars <- c(A = .8, B = -.3, sigma_A = -1, sigma_B = 1)
#' myfn(pars = pars[c(1, 3)], fixed = pars[c(2, 4)])
#' 
#' # Assume same sigma parameter for both A and B
#' # sigma is assumed to be passed on log scale
#' mu <- c(A = 0, B = 0)
#' myfn <- constraintL2(mu, sigma = "sigma")
#' pars <- c(A = .8, B = -.3, sigma = 0)
#' myfn(pars = pars)
#' 
#' @export
constraintL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {

  
  # Aktuell zu kompliziert aufgesetzt. Man sollte immer die komplette Hessematrix/Gradient
  # auswerten und dann die Elemente streichen, die in fixed sind!
  
  
  estimateSigma <- ifelse(is.character(sigma), TRUE, FALSE)
  if (length(sigma) > 1 & length(sigma) < length(mu))
    stop("sigma must either have length 1 or at least length equal to length of mu.")
  
  ## Augment sigma if length = 1
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  if (is.null(names(sigma)))
    names(sigma) <- names(mu)
  if (!is.null(names(sigma)) & !all(names(mu) %in% names(sigma)))
    stop("Names of sigma and names of mu do not match.")
  
  ## Bring sigma in correct order (no matter if character or numeric)
  sigma <- sigma[names(mu)]
  
  controls <- list(mu = mu, sigma = sigma, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Import from controls
    mu <- controls$mu
    sigma <- controls$sigma
    attr.name <- controls$attr.name
    nmu <- length(mu)
    
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
      
      
      pars <- c(p, fixed)[names(mu)]
      p1 <- setdiff(intersect(names(mu), names(p)), names(fixed))
      
      # if estimate sigma, produce numeric sigma vector from the parameters provided in p and fixed
      if (estimateSigma) {
        sigmapars <- sigma
        sigma <- exp(c(p, fixed)[sigma])
        names(sigma) <- names(mu)
        Jsigma <- do.call(cbind, lapply(unique(sigmapars), function(s) {
          (sigmapars == s)*sigma
        }))
        colnames(Jsigma) <- unique(sigmapars)
        rownames(Jsigma) <- names(sigma)
        p2 <- setdiff(intersect(unique(sigmapars), names(p)), names(fixed))
      }
      
      # Compute constraint value and derivatives
      val <- sum((pars - mu)^2/sigma^2) + estimateSigma * sum(log(sigma^2))
      val.p <- 2*(pars - mu)/sigma^2
      val.sigma <- -2*(pars-mu)^2/sigma^3 + 2/sigma
      val.p.p <- diag(2/sigma^2, nmu, nmu); colnames(val.p.p) <- rownames(val.p.p) <- names(mu)
      val.p.sigma <- diag(-4*(pars-mu)/sigma^3, nmu, nmu); colnames(val.p.sigma) <- rownames(val.p.sigma) <- names(mu)
      val.sigma.sigma <- diag(6*(pars-mu)^2/sigma^4 - 2/sigma^2, nmu, nmu); colnames(val.sigma.sigma) <- rownames(val.sigma.sigma) <- names(mu)
      
      # Multiply with Jacobian of sigma vector if estimate sigma
      if (estimateSigma) {
        val.sigma.sigma <- t(Jsigma) %*% val.sigma.sigma %*% Jsigma + diag((t(val.sigma) %*% Jsigma)[1,], ncol(Jsigma), ncol(Jsigma))
        val.sigma <- (val.sigma %*% Jsigma)[1,]
        val.p.sigma <- (val.p.sigma %*% Jsigma)
      }
      
      
      gr <- hs <- NULL
      if (deriv) {
        # Produce output gradient and hessian
        gr <- rep(0, length(p)); names(gr) <- names(p)
        hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
        
        # Set values in gradient and hessian
        gr[p1] <- val.p[p1]
        hs[p1, p1] <- val.p.p[p1, p1]
        if (estimateSigma) {
          gr[p2] <- val.sigma[p2]
          hs[p1, p2] <- val.p.sigma[p1, p2]
          hs[p2, p1] <- t(val.p.sigma)[p2, p1]
          hs[p2, p2] <- val.sigma.sigma[p2, p2]
        }
        
        # Multiply with derivatives of incoming parameter
        dP <- attr(p, "deriv")
        if (!is.null(dP)) {
          gr <- as.vector(gr %*% dP); names(gr) <- colnames(dP)
          hs <- t(dP) %*% hs %*% dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
        }
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
    if (!is.null(conditions) && !condition %in% conditions) 
      return()
    if (is.null(conditions) && !condition %in% names(prediction))
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")
    
    
    
    # Divide parameter into data point and rest
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    
    
    # Get predictions and derivatives at time point
    time.index <- which(prediction[[condition]][,"time"] == time)
    if (length(time.index) == 0) 
      stop("datapointL2() requests time point for which no prediction is available. Please add missing time point by the times argument in normL2()")
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
    attr(out, "prediction") <- pred
    
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
      
      gr <- hs <- NULL
      if (deriv) {
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
  
  # Extract BLOQ part from nout
  is.bloq <- nout$bloq
  nout.bloq <- nout[is.bloq, , drop = FALSE]
  
  # Drop BLOQ part from nout
  nout <- nout[!is.bloq, , drop = FALSE]
  obj <- sum(nout$weighted.residual^2) + sum(-2*log(pnorm(-nout.bloq$weighted.residual)))
  
  grad <- NULL
  hessian <- NULL
  derivs <- attr(nout, "deriv")
  if (!is.null(derivs)) {
    
    # Extract BLOQ part from derivs
    derivs.bloq <- derivs[is.bloq, , drop = FALSE]
    # Drop BLOQ part from derivs
    derivs <- derivs[!is.bloq, , drop = FALSE]

    if (nrow(derivs) > 0) {
      
      nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
      sens <- as.matrix(derivs[, -(1:2), drop = FALSE])
      
      res <- nout$residual
      sigma <- nout$sigma
      
      grad <- as.vector(2*matrix(res/sigma^2, nrow = 1) %*% sens)
      names(grad) <- colnames(sens)
      hessian <- 2*t(sens/sigma) %*% (sens/sigma)
      
    }
    
    if (nrow(derivs.bloq) > 0) {
      
      nout.bloq[is.na(nout.bloq$sigma)] <- 1
      sens.bloq <- as.matrix(derivs.bloq[, -(1:2), drop = FALSE])
      
      Phi <- pnorm(-nout.bloq$weighted.residual)
      G <- dnorm(-nout.bloq$weighted.residual)
      res <- nout.bloq$residual
      sigma <- nout.bloq$sigma
      
      
      grad.bloq <- as.vector(matrix(2*G/(Phi*sigma), nrow = 1) %*% sens.bloq)
      names(grad.bloq) <- colnames(sens.bloq)
      
      X1 <- sens.bloq*(G/(Phi*sigma))^2
      X2 <- sens.bloq*(res*G/(Phi*sigma^3))
      hessian.bloq <- 2 * t(X1) %*% sens.bloq - 2 * t(X2) %*% sens.bloq
      
      if (is.null(grad) & is.null(hessian)) {
        grad <- grad.bloq
        hessian <- hessian.bloq
      } else {
        grad <- grad + grad.bloq
        hessian <- hessian + hessian.bloq
      }
      
    }
    
    
    
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
#' @importFrom stats pnorm dnorm
nll <- function(nout) {
  
  # Extract BLOQ part from nout
  is.bloq <- nout$bloq
  nout.bloq <- nout[is.bloq, , drop = FALSE]
  
  # Drop BLOQ part from nout
  nout <- nout[!is.bloq, , drop = FALSE]
  
  obj <- sum(nout$weighted.residual^2) + sum(log(2*pi*nout$sigma^2)) + 
    sum(-2*log(pnorm(-nout.bloq$weighted.residual)))
  grad <- NULL
  hessian <- NULL
  
  
  
  derivs <- attr(nout, "deriv")
  derivs.err <- attr(nout, "deriv.err")
  if (!is.null(derivs) & !is.null(derivs.err)) {
    
    # Extract BLOQ part from derivs
    derivs.bloq <- derivs[is.bloq, , drop = FALSE]
    derivs.err.bloq <- derivs.err[is.bloq, , drop = FALSE]
    # Drop BLOQ part from derivs
    derivs <- derivs[!is.bloq, , drop = FALSE]
    derivs.err <- derivs.err[!is.bloq, , drop = FALSE]
    
    if (nrow(derivs) > 0 & nrow(derivs.err) > 0) {
    
      # Get sensitivities: sens = dres/dp, sens.err = dsigma/dp
      sens <- as.matrix(derivs[, -(1:2), drop = FALSE])
      sens.err <- as.matrix(derivs.err[, -(1:2), drop = FALSE])
      
      res <- nout$residual
      sigma <- nout$sigma
      
      # Compute gradient
      grad <- as.vector(2*matrix(res/sigma^2, nrow = 1) %*% sens -
                          2*matrix(res^2/sigma^3, nrow = 1) %*% sens.err +
                          2*matrix(1/sigma, nrow = 1) %*% sens.err)
      names(grad) <- colnames(sens)
      
      # Compute hessian
      X1 <- (1/sigma)*sens - (res/sigma^2)*sens.err
      X2 <- (res/sigma^2)*sens.err
      X3 <- (1/sigma)*sens.err
      
      hessian <- 2 * t(X1) %*% X1 #+ 4 * t(X2) %*% X2 - 2 * t(X3) %*% X3
      
    }
    
    if (nrow(derivs.bloq) > 0 & nrow(derivs.err.bloq) > 0) {
      
      # Get sensitivities: sens = dres/dp, sens.err = dsigma/dp
      sens.bloq <- as.matrix(derivs.bloq[, -(1:2), drop = FALSE])
      sens.err.bloq <- as.matrix(derivs.err.bloq[, -(1:2), drop = FALSE])
      
      Phi <- pnorm(-nout.bloq$weighted.residual)
      G <- dnorm(-nout.bloq$weighted.residual)
      res <- nout.bloq$residual
      sigma <- nout.bloq$sigma
      
      # Compute gradient
      grad.bloq <- as.vector(matrix(2*G/(Phi*sigma), nrow = 1) %*% sens.bloq) -
        as.vector(matrix(2*G*res/(Phi*sigma^2), nrow = 1) %*% sens.err.bloq)
      names(grad.bloq) <- colnames(sens.bloq)
      
      # Compute hessian
      X <- -(1/sigma)*sens.bloq + (res/sigma^2)*sens.err.bloq
      X1 <- (2*G^2/Phi^2)*X
      X2 <- (2*G*res/(Phi*sigma))*X
      
      hessian.bloq <- t(X1) %*% X  - t(X2) %*% X
      
      if (is.null(grad) & is.null(hessian)) {
        grad <- grad.bloq
        hessian <- hessian.bloq
      } else {
        grad <- grad + grad.bloq
        hessian <- hessian + hessian.bloq
      }
      
    }
    
    
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

