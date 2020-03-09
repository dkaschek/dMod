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
#' 
#' 
#' @param forcings TO BE FILLED BY DANIEL K
#' @param iiv Example: c("ka", "ETA_EMAX"). \cr 
#'   Vector with names which are individualized per condition
#' @param conditional Example: data.frame(parname = "GR", covname = "SEX", covvalue = "1", stringsAsFactors = FALSE).\cr
#'   * covname can relate to any parameter in the condition.grid of the data. \cr
#'   * covvalue is the value of this variable to use for individualization
#'
#' @param fixed.grid data.frame(parname, partask, ids...) Lookup table for fixed parameters
#' @param nauxtimes additional simulation times
#' @param cores to parallelize over conditions not over fits
#' 
#' 
#' @return Object of class \code{obsfn}, i.e. a function 
#' \code{obj(..., fixed, deriv, conditions, env)} that returns an objective list,
#' \link{objlist}.
#' @details Objective functions can be combined by the "+" operator, see \link{sumobjfn}.
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL, attr.name = "data") {

  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times)) timesD <- sort(union(times, timesD))

  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  if (!all(data.conditions %in% x.conditions)) 
    stop("The prediction function does not provide predictions for all conditions in the data.")
  e.conditions <- names(attr(errmodel, "mappings"))
  
  controls <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions, data.conditions))

  # might be necessary to "store" errmodel in the objective function (-> runbg)
  force(errmodel)  
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = controls$conditions, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Generate output template
    pars_out <- colnames(getDerivs(as.parvec(pouter)))
   
    # Import from controls
    timesD <- controls$times
    attr.name <- controls$attr.name
    
    # Create new environment if necessary
    if (is.null(env)) env <- new.env()
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) {
      err <- NULL
      if ((!is.null(errmodel) & is.null(e.conditions)) | (!is.null(e.conditions) && (cn %in% e.conditions))) 
        err <- errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
      nll(res(data[[cn]], prediction[[cn]], err[[cn]]), pars = pouter, deriv = deriv)
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


#' Compute the negative log-likelihood
#' 
#' Gaussian Log-likelihood. Supports NONMEM-like BLOQ handling methods M1, M3 and M4 and estimation of error models
#' 
#' @param nout data.frame (result of [res]) or object of class [res].
#' @param pars Example named vector of outer parameters to construct the objlist
#' @param deriv TRUE or FALSE
#' @param opt.BLOQ see [normIndiv]
#' @param opt.hessian see [normIndiv]
#' 
#' @md
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
#' @export
nll <- function(nout, pars, deriv, opt.BLOQ = "M3", opt.hessian = c(
  ALOQ_part1 = TRUE, ALOQ_part2 = TRUE, ALOQ_part3 = TRUE,
  BLOQ_part1 = TRUE, BLOQ_part2 = TRUE, BLOQ_part3 = TRUE,
  PD = TRUE  # enforce Hessian to be positive semidefinite, by setting nearest negative eigenvalues to zero
  )) {
  
  # Split residuals into ALOQ and BLOQ
  is.bloq   <- nout$bloq
  nout.bloq <- nout[is.bloq, , drop = FALSE]
  nout.aloq <- nout[!is.bloq, , drop = FALSE]
  
  # Handle derivs
  derivs          <- attr(nout, "deriv")
  derivs.bloq     <- derivs[is.bloq, , drop = FALSE]
  derivs.aloq     <- derivs[!is.bloq, , drop = FALSE]
  derivs.err      <- attr(nout, "deriv.err")
  derivs.err.bloq <- derivs.err[is.bloq, , drop = FALSE]
  derivs.err.aloq <- derivs.err[!is.bloq, , drop = FALSE]
  
  # Apply nll
  mywrss <- init_empty_objlist(pars, deriv = deriv)
  nll_ALOQ_result <- NULL
  if (!all(is.bloq))
    nll_ALOQ_result <- nll_ALOQ(nout.aloq, derivs.aloq, derivs.err.aloq, opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian)
  mywrss <- mywrss + nll_ALOQ_result
  if (any(is.bloq) && (!opt.BLOQ == "M1"))
    mywrss <- mywrss + nll_BLOQ(nout.bloq, derivs.bloq, derivs.err.bloq, opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian)
  
  attr(mywrss, "chisquare") <- attr(nll_ALOQ_result, "chisquare")
  
  mywrss
}



#' Non-linear log likelihood for the ALOQ part of the data
#' @md
#' @param nout output of [res()]
#' @param derivs,derivs.err attributes of output of [res()]
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector to include or exclude
#'   various non-convex summands of the hessian matrix
#' @importFrom stats pnorm dnorm
nll_ALOQ <- function(nout,
                     derivs,
                     derivs.err,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(
                       # ALOQ (Above limit of quantification): Three parts of the hessian which can be obtained with first order derivs
                       ALOQ_part1 = TRUE,
                       ALOQ_part2 = TRUE,
                       ALOQ_part3 = TRUE
                     )) {
  
  # opt.BLOQ: only M4BEAL changes something, all others are implicit, since they do not change anything in the ALOQ part
  # opt.hessian ALOQ_part1-ALOQ_part3
  
  # .. Residual terms ----#
  wr <- nout$weighted.residual
  w0 <- nout$weighted.0
  s  <- nout$sigma
  
  # .. Compute objective value ----#
  chisquare <- obj <- sum(wr^2)
  obj <- obj + sum(log(2*pi*s^2))
  if (opt.BLOQ %in% "M4BEAL")
    obj <- obj + 2 * sum(stats::pnorm(w0, log.p = TRUE))
  
  grad <- NULL
  hessian <- NULL
  if (!is.null(derivs) && nrow(derivs) > 0) {
    # .. Sensitivities terms ----#
    #   sens = dres/dp = dx/dp, sens.err = dsigma/dp
    dxdp <- as.matrix(derivs[, -(1:2), drop = FALSE])
    dsdp <- 0 * dxdp
    if (!is.null(derivs.err))
      dsdp <- as.matrix(derivs.err[, -(1:2), drop = FALSE])
    dwrdp <- 1/s*dxdp - wr/s*dsdp
    dw0dp <- 1/s*dxdp - w0/s*dsdp
    dlogsdp <- (1/s)*dsdp # dlogsig.dp
    G_by_Phi <- function(w) exp(stats::dnorm(w, log = TRUE)- stats::pnorm(w, log.p = TRUE))
    # .. 2nd Sensitivity terms ----#
    #   interaction terms of second derivative. d2adb2 means second derivative: d^2a/db^2
    #   d2wrdp2 does not solely consist of second derivatives but via the prefactors also of some combinations of first order derivs.
    # * These are the equations, but since we need to sum over residuals, they will be inserted directly
    # d2wrdp2 <- - 1 / (s^2) * (dxdp * dsdp + dsdp * dxdp) +
    #               2 * wr/(s^2) * dsdp * dsdp # - wr/s* d2sdp2 + 1/s * d2xdp2
    
    # .. Compute gradient ----#
    grad <- as.vector(2*matrix(wr, nrow = 1) %*% dwrdp + 2*apply(dlogsdp,2, sum))
    if (opt.BLOQ %in% "M4BEAL")
      grad <- grad + as.vector((2 * matrix(G_by_Phi(w0), nrow = 1)) %*% (dw0dp))
    names(grad) <- colnames(dxdp)
    
    # .. Compute hessian ----#
    # >>>> All equations were double-checked, they should be fine. (Dont touch or read them, D2!) <<<<<<<<<<<
    hessian <- matrix(0, nrow = ncol(dwrdp), ncol = ncol(dwrdp), dimnames = list(colnames(dwrdp), colnames(dwrdp)))
    hessian <- hessian + 2 * t(dwrdp) %*% dwrdp # - 2 * t(dsig.dp) %*% dsig.dp # + 2. sens
    
    if (opt.hessian["ALOQ_part1"])  # "interaction"-terms of d2wrdp2
      hessian <- hessian + 2 * (t(-wr/s^2 * dxdp) %*% dsdp + t(-wr/s^2 * dsdp) %*% dxdp) #slightly concave
    if (opt.hessian["ALOQ_part2"])  # "interaction"-terms of d2wrdp2
      hessian <- hessian + 2 * t(2 * wr^2/(s^2) * dsdp)%*%dsdp # + 2nd order derivs of logsig and wr
    if (opt.hessian["ALOQ_part3"]) # The non-convex contribution by log(sigma)
      hessian <- hessian - 2 * t(dlogsdp) %*% dlogsdp
    
    
    if (opt.BLOQ %in% "M4BEAL") {
      hessian <- hessian + 2 * t((-w0 * G_by_Phi(w0) - G_by_Phi(w0)^2) * dw0dp) %*% dw0dp # 1st order terms
      hessian <- hessian + 2 * t(G_by_Phi(w0) * (-1)/(s^2) * dxdp ) %*% dsdp + 2 * t(G_by_Phi(w0) * (-1)/(s^2) * dsdp ) %*% dxdp # 2nd order terms
      if (opt.hessian["ALOQ_part1"]) # The other term from second sensitivities of wr
        hessian <- hessian + 2 * t(2 * G_by_Phi(w0) * w0/(s^2) * dsdp)%*%dsdp
    }
  }
  
  out <- objlist(value = obj, gradient = grad, hessian = hessian)
  attr(out, "chisquare") <- chisquare
  out
}



#' Non-linear log likelihood for the BLOQ part of the data
#' @md
#' @param nout.bloq The bloq output of [res()]
#' @param derivs.bloq,derivs.err.bloq attributes of output of [res()]
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector to include or exclude
#'   various summands of the hessian matrix
#' @importFrom stats pnorm dnorm
nll_BLOQ <- function(nout.bloq,
                     derivs.bloq,
                     derivs.err.bloq,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(
                       BLOQ_part1 = TRUE,
                       BLOQ_part2 = TRUE,
                       BLOQ_part3 = TRUE
                     )) {
  
  # .. Checks -----
  if (opt.BLOQ %in% c("M4NM", "M4BEAL") & any(nout.bloq$value < 0))
    stop("M4-Method cannot handle LLOQ < 0. Possible solutions:
      * Use M3 which allows negative LLOQ (recommended)
      * If you are working with log-transformed DV, exponentiate DV and LLOQ\n")
  
  # .. Residuals and sensitivities ----#
  wr <- nout.bloq$weighted.residual
  w0 <- nout.bloq$weighted.0
  s  <- nout.bloq$sigma
  
  # .. Compute objective value ----#
  if (opt.BLOQ == "M3"){
    objvals.bloq <- -2*stats::pnorm(-wr, log.p = TRUE)
  }
  if (opt.BLOQ %in% c("M4NM", "M4BEAL")){
    objvals.bloq <- -2*log(1 - stats::pnorm(wr) / stats::pnorm(w0))
    # .... catch numerically problematic cases ----#
    # The problematic region can be approximated by a parabola, intercept and linear coefficient depend on LOQ/s
    intercept = ifelse(log(w0-wr) > 0, 1.8, -1.9 * log(w0-wr) +0.9)
    lin = ifelse(log(w0-wr) > 0, 0.9, 0.5 )
    objvals.bloq[!is.finite(objvals.bloq)] <-  (intercept + lin * w0 + 0.95 * w0^2)[!is.finite(objvals.bloq)]
  }
  
  obj.bloq <- sum(objvals.bloq)
  grad.bloq <- NULL
  hessian.bloq <- NULL
  
  if (!is.null(derivs.bloq) && nrow(derivs.bloq) > 0){
    # .. Sensitivities ----#
    #   sens = dres/dp = dx/dp, sens.err = dsigma/dp = dsdp
    dxdp <- as.matrix(derivs.bloq[, -(1:2), drop = FALSE])
    dsdp <- 0 * dxdp
    if (!is.null(derivs.err.bloq))
      dsdp <- as.matrix(derivs.err.bloq[, -(1:2), drop = FALSE])
    dwrdp <- 1/s*dxdp - wr/s*dsdp
    dw0dp <- 1/s*dxdp - w0/s*dsdp
    dlogsdp <- (1/s)*dsdp # dlogsig.dp
    G_by_Phi <- function(w1, w2 = w1) exp(stats::dnorm(w1, log = TRUE) - stats::pnorm(w2, log.p = TRUE))
    
    # .. 2nd Sensitivities ----#
    #   interaction terms of second derivative. d2adb2 means second derivative: d^2a/db^2
    #   d2wrdp2 does not solely consist of second derivatives but via the prefactors also of some combinations of first order derivs.
    # * These are the equations, but since we need to sum over residuals, they will be inserted directly
    # d2wrdp2 <- - 1 / (s^2) * (dxdp * dsdp + dsdp * dxdp) +
    #               2 * wr/(s^2) * dsdp * dsdp # - wr/s* d2sdp2 + 1/s * d2xdp2
    
    # .. Compute gradient ----#
    if (opt.BLOQ == "M3"){
      grad.bloq <- -2 * as.vector(matrix( G_by_Phi(-wr), nrow = 1) %*% (-dwrdp)) # minus sign of -2logLikelihood and -dwrdp cancel # [] clean this formula
    }
    if (opt.BLOQ %in% c("M4NM", "M4BEAL")){
      grad.bloq <-             as.vector(matrix(2 / (1/G_by_Phi(wr,w0) - 1/G_by_Phi(wr,wr)), nrow = 1) %*% dwrdp)
      grad.bloq <- grad.bloq - as.vector(matrix(2 / (1/G_by_Phi(w0,w0) - 1/G_by_Phi(w0,wr)), nrow = 1) %*% dw0dp)
      grad.bloq <- grad.bloq + as.vector(matrix(2 * G_by_Phi(w0), nrow = 1) %*% dw0dp)
    }
    names(grad.bloq) <- colnames(dxdp)
    
    # .. Compute hessian ----
    if (opt.BLOQ %in% "M3") {
      hessian.bloq <- matrix(0, nrow = ncol(dxdp), ncol = ncol(dxdp), dimnames = list(colnames(dxdp), colnames(dxdp)))
      if (opt.hessian["BLOQ_part1"])
        hessian.bloq <- hessian.bloq + 2 * t((-wr * G_by_Phi(-wr) + G_by_Phi(-wr)^2) * dwrdp) %*% dwrdp # 1st order terms
      if (opt.hessian["BLOQ_part2"]){
        hessian.bloq <- hessian.bloq - 2 * t(G_by_Phi(-wr) * (+1)/(s^2) * dxdp) %*% dsdp # 2nd order terms (+1) because the original term of d2wrdp2 contains -1/s^2, so (d2-wrdp) contains +1/s^2
        hessian.bloq <- hessian.bloq - 2 * t(G_by_Phi(-wr) * (+1)/(s^2) * dsdp) %*% dxdp
      }
      if (opt.hessian["BLOQ_part3"])
        hessian.bloq <- hessian.bloq - 2 *  t(G_by_Phi(-wr) * (2 * (-wr))/(s^2) * dsdp)%*%dsdp
    }
    
    # .. M4 hessian ----
    if (opt.BLOQ %in% c("M4NM", "M4BEAL")) {
      d_dp_sq <- function(A, w = wr, sign = 1) {
        # @details This function performs the multiplication of A with quadratic first order derivs: A*t(dw/dp)%*%(dw/dp)
        # @param A: sth to multiply the derivs with
        # @param w weighted residual: practically choice between wr or w0
        # @param sign: +1 or -1, e.g. in H(-2*LLBM3), wr appears as -wr, therefore, the minus signs have to be propagated correctly. (In A, however, everything has to specified by hand).
        dwdp <- 1/s*dxdp - w/s*dsdp
        out <- matrix(0, nrow = ncol(dxdp), ncol = ncol(dxdp), dimnames = list(colnames(dxdp), colnames(dxdp)))
        out <- out + t(A * dwdp) %*% dwdp
      }
      d2_dp2 <- function(A, w = wr, sign = 1) {
        # @details This function performs the multiplication of A with the second derivative terms: A*d^2w/dp^2
        # @inheritParams d_dp_sq
        out <- matrix(0, nrow = ncol(dxdp), ncol = ncol(dxdp), dimnames = list(colnames(dxdp), colnames(dxdp)))
        out <- out + t(A * (-1*sign)/(s^2) * dxdp) %*% dsdp
        out <- out + t(A * (-1*sign)/(s^2) * dsdp) %*% dxdp
        out <- out + t(A * (2 * (w*sign))/(s^2) * dsdp) %*% dsdp
      }
      # In hessian(BM4), one often gets an expression G(w0 or wr)/(Phi(w0)-Phi(wr))
      # This function is replaces numerically problematic regions by an approximation
      # @params wn,w0,wr "Weighted residual Numerator/Denominator 1/2"
      stable <- function(wn, w0, wr) {
        
        if(!(identical(wn, w0) | identical(wn, wr)))
          stop("The first argument wn needs to be identical to either the second or third")
        
        out <- stats::dnorm(wn)/(stats::pnorm(w0)-stats::pnorm(wr))
        # two possible cases with different asymptotics: wn == wd1 or wn == wd2
        if (identical(wn, w0)){
          out[is.infinite(out)] <- 0
          return(out)
        }
        if (identical(wn, wr)){
          out[is.infinite(out)] <- 1/(w0-wr) + wr # This formula was found out "by hand", I didn't really search for an analytic justification. If you want, you can try l'Hospitalizing it.
          return(out)
        }
      }
      
      part1 <- d_dp_sq(-wr * stable(wr,w0,wr)) +
        d2_dp2(stable(wr,w0,wr)) -
        (d_dp_sq(-w0 * stable(w0,w0,wr), w = w0) +
           d2_dp2(stable(w0,w0,wr), w = w0))
      part1 <- 2 * part1
      
      part2 <- stable(wr,w0,wr) * dwrdp - stable(w0,w0,wr) * dw0dp
      part2 <- -2 * t(part2) %*% part2
      
      part3 <- d_dp_sq(-w0 * G_by_Phi(w0) - (G_by_Phi(w0))^2, w = w0) + d2_dp2(G_by_Phi(w0), w = w0)
      part3 <- 2 * part3
      
      hessian.bloq <- matrix(0, nrow = ncol(dxdp), ncol = ncol(dxdp), dimnames = list(colnames(dxdp), colnames(dxdp)))
      if (opt.hessian["BLOQ_part1"])
        hessian.bloq <- hessian.bloq + part1
      if (opt.hessian["BLOQ_part2"])
        hessian.bloq <- hessian.bloq + part2
      if (opt.hessian["BLOQ_part3"])
        hessian.bloq <- hessian.bloq + part3
    }
  }
  
  out <- objlist(value = obj.bloq, gradient = grad.bloq, hessian = hessian.bloq)
  return(out)
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

  what <- intersect(c("value", "gradient", "hessian"), c(names(out1), names(out2)))
  
  add_vector <- function(a,b) {
    # add vector b to a by names
    i <- intersect(names(a), names(b))
    a[i] <- a[i] + b[i]
    a}
  add_matrix <- function(a,b) {
    i <- intersect(rownames(a), rownames(b))
    a[i,i] <- a[i,i] + b[i,i]
    a}
  
  gn1 <- names(out1$gradient)
  gn2 <- names(out2$gradient)
  
  one_includes_two <- all(gn2 %in% gn1) 
  two_includes_one <- all(gn1 %in% gn2)
  neither_included <- !(one_includes_two | two_includes_one)
  
  out12 <- lapply(what, function(w) {
    v1 <- out1[[w]]
    v2 <- out2[[w]]
    if (w == "value") 
      return(v1 + v2)
    if (w == "gradient"){
      if (neither_included) return(add_vector(add_vector(setNames(rep(0, length(union(gn1, gn2))), union(gn1, gn2)),v1),v2))
      if (one_includes_two) return(add_vector(v1,v2))
      if (two_includes_one) return(add_vector(v2,v1))
    }
    if (w == "hessian") {
      if (neither_included) return(add_matrix(add_matrix(matrix(0, length(union(gn1,gn2)),length(union(gn1,gn2)),
                                                       dimnames = list(union(gn1,gn2), union(gn1,gn2))
                                                       ),v1),v2))
      if (one_includes_two) return(add_matrix(v1,v2))
      if (two_includes_one) return(add_matrix(v2,v1))
    }
  })
  names(out12) <- what
  
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

