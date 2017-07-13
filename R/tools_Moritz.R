#######################################################################################################################################################################

parametersToPoints <- function(startpars, n, timepoints, sigma, truepars, par_names, model, choice = 1:length(timepoints[[1]])){
  #variable to collect all parameters with predefined size
  par_total <- list(startpars)
  #set data points which will replace the initial parmeters
  datapoints_unsort <- wide2long((model)(times = timepoints[[1]], truepars))
  datapoints <- do.call(rbind,lapply(1:length(timepoints[[1]]), function(i) datapoints_unsort[which(datapoints_unsort$time==timepoints[[1]][i])[1],] ))
  datapoints <- rbind(datapoints,do.call(rbind,lapply(1:length(timepoints[[1]]), function(i) datapoints_unsort[which(datapoints_unsort$time==timepoints[[1]][i])[2],] )))
  
  #sort out unwanted points like the event time
  datapoints_obs1 <- datapoints[1:(dim(datapoints)[1]/2),]
  datapoints_obs1 <- do.call(rbind,lapply(1:length(timepoints[[1]]), function(i){if(timepoints[[2]][[i]]=="obs1"){return(datapoints_obs1[i,])}}))
  
  datapoints_obs2 <- datapoints[((dim(datapoints)[1]/2)+1):(dim(datapoints)[1]),]
  datapoints_obs2 <- do.call(rbind,lapply(1:length(timepoints[[1]]), function(i){if(timepoints[[2]][[i]]=="obs2"){return(datapoints_obs2[i,])}}))
  
  datapoints <- rbind(datapoints_obs1,datapoints_obs2)
  datapoints <- do.call(rbind, lapply(1:dim(datapoints)[1], function(i) datapoints[choice[i],] ))
  
  datapoints_choice <- NULL
  for(j in 1:length(choice)){#j = number of points
    #error output
    if( !(length(timepoints[[1]])==length(choice)) ){
      return("Choice should have the same length as timepoints!")
      break}
    
    datapoints_choice$sigma     <- NULL
    datapoints_choice           <- rbind(datapoints_choice,datapoints[j,])
    datapoints_choice$sigma     <- sigma
    datapoints_choice_list      <- as.datalist(datapoints_choice)
    
    obj <- normL2(datapoints_choice_list, model)
    
    #fit and overwrite parameters
    par_total[[j+1]] <- do.call(rbind,
                                mclapply(1:n,
                                         function(i){
                                           myfit<-trust(obj, parinit = structure(as.numeric(par_total[[j]][i,]), names = par_names), rinit = 1, rmax = 10)
                                           if(myfit$converged && (myfit$value<=1/(10^5*sigma)) ){return(myfit$argument)}
                                           else{return(structure(c(0,0,0,0), names = par_names))}
                                         },
                                         mc.cores = 24,
                                         mc.preschedule = FALSE ))
    
  }
  par_total[[length(choice)+2]]=datapoints
  return(par_total)
}

#######################################################################################################################################################################

parametersToPointsVar <- function(startpars, sigma, truepars, par_names, model, timeframe = c(0,50,.1), tol = 1/sigma[1], cores = 24){
  n <- dim(startpars)[1]
  
  #for-loop setup
  par_total     <- list(startpars)
  par_total_del <- list(startpars)
  delold        <- c(n+1)
  timepoints    <- list(c(),c())
  timepoint     <- list(c(),c())
  derivs_save   <- list()
  
  for(j in 1:length(truepars)){#j = number of points
    
    #find sensitivity maximum
    derivs <- list()
    times <- seq(from = timeframe[1], to = timeframe[2], by = timeframe[3])
    par              <- par_total_del[[j]]
    meanvals         <- Reduce("+", mclapply(1:dim(par_total_del[[j]])[1], function(i) (wide2long((model)(times, par[i,]))$value), mc.cores = cores, mc.preschedule = TRUE ))/dim(par_total_del[[j]])[1]
    derivs[[1]]      <- Reduce("+", mclapply(1:dim(par_total_del[[j]])[1], function(i) (meanvals-wide2long((model)(times, par[i,]))$value)^2, mc.cores = cores, mc.preschedule = TRUE ))/(dim(par_total_del[[j]])[1]-1) #with mean value
    derivs[[2]]      <- rep(times,2)
    
    
    #take the parameter with the highest sensitivity first, then the one with the second highest and s.o .... 
    sigma_ratio <- sigma[1]/sigma[2]
    
    observables <- names(attr(attr(model, "mappings")[[1]], "equations"))
    observables_max <- do.call(rbind,lapply(1:length(observables), function(i) max(derivs[[1]][((i-1)*length(times)+1):(i*length(times))]) ))
    timepoint[[1]] <- derivs[[2]][which.max(derivs[[1]])]
    timepoint[[2]] <- observables[which.max(observables_max)]
    
    timepoints[[1]][j] <- timepoint[[1]]
    timepoints[[2]][j] <- timepoint[[2]]
    
    #set data points which will replace the initial parmeters
    datapoints_unsort <- wide2long((model)(times = timepoints[[1]], truepars))
    datapoints <- do.call(rbind,lapply(1:length(timepoints[[1]]), function(i) datapoints_unsort[which(datapoints_unsort$time==timepoints[[1]][i])[1],] )) #delete unwanted t=30 for obs1
    datapoints <- rbind(datapoints,do.call(rbind,lapply(1:length(timepoints[[1]]), function(i) datapoints_unsort[which(datapoints_unsort$time==timepoints[[1]][i])[2],] ))) #delete unwanted t=30 for obs2
    
    #sort out unwanted points
    datapoints_obs <- data.frame(time = c(), name = factor(levels = observables), value = c(),condition = factor(levels = "control"))
    for(k in 1:length(observables)){
      datapoints_dummy <- datapoints[((k-1)*(dim(datapoints)[1]/length(observables))+1):(k*(dim(datapoints)[1]/length(observables))),]
      datapoints_obs <- rbind(datapoints_obs,do.call(rbind,lapply(1:length(timepoints[[1]]), function(i){if(timepoints[[2]][[i]]==observables[k]){return(datapoints_dummy[i,])}})))
    }
    
    datapoints <- datapoints_obs
    datapoints$sigma <- NULL
    
    datapoints$sigma <- sigma[1]
    datapoints_list  <- as.datalist(datapoints)
    
    obj <- normL2(datapoints_list, model)
    
    #fit and overwrite parameters
    par_total[[j+1]] <- do.call(rbind,
                                mclapply(1:dim(par_total[[j]])[1],
                                         function(i){
                                           myfit <- trust(obj, parinit = structure(as.numeric(par_total[[j]][i,]), names = par_names), rinit = .1, rmax = 10)
                                           if(myfit$converged && (myfit$value<=tol/(sigma[1])) ){return(myfit$argument)}
                                           else{return(structure(rep(0, length(par_names)), names = par_names))}
                                         },
                                         mc.cores = cores,
                                         mc.preschedule = FALSE ))
    
    #NaNs need to be removed
    del <- intersect(which(as.data.frame(par_total[[j+1]])$import_Tca == 0), which(as.data.frame(par_total[[j+1]])$export_Tca_baso == 0))
    del <- union(union(del,c(n+1)),delold) #n+1 is implemented because del must not be empty
    par_total_del[[j+1]] <- as.data.frame(par_total[[j+1]])
    
    par_total_del[[j+1]] <- par_total_del[[j+1]][-del,]
    delold <- del
    
    print(par_total_del)
  }
  
  par_total_del[[length(truepars)+2]] <- datapoints
  par_total_del[[length(truepars)+3]] <- timepoints
  par_total_del[[length(truepars)+4]] <- derivs_save
  return(par_total_del)
}

#######################################################################################################################################################################

norm <- function(x) sqrt(sum(x^2))
###
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
###
Pdata <- function(times, obs_names, par, model, thresh = 1e-3){
  sensnames      <- do.call(rbind, lapply(1:length(times), function(i) paste0(rep(obs_names[i],length(par)), rep(".",length(par)), names(par) ))) #parinit must have names
  sensnames      <- cbind(sensnames,rep("time",length(par)))
  
  sens           <- do.call(rbind, lapply(1:length(times), function(i) getDerivs((model)(times[i],par))[[1]][,unlist(sensnames[i,])]))
  sens           <- do.call(rbind, lapply(1:length(times), function(i) sens[which(sens[,"time"]==times[i]),-(length(times)+1)] )) #sort to correct order again, deletes unwanted time=30
  
  sens <- sens
  
  colnames(sens) <- c(names(par))
  rownames(sens) <- as.character(do.call(cbind, lapply(1:length(times), function(i) paste0("d_",i) )))
  
  mysvd <- svd(sens)
  mydiag <- (mysvd$d>=thresh)*mysvd$d
  
  if(sum(mydiag == 0)>=1){print(paste("sensitivity suppressed:","d=",mysvd$d[1],mysvd$d[2],mysvd$d[3],mysvd$d[4]))}
  
  sens_inv <- pseudoinverse(sens,thresh)
  
  colnames(sens_inv) <- as.character(do.call(cbind, lapply(1:length(times), function(i) paste0("d_",i) )))
  rownames(sens_inv) <- c(names(par))
  
  datapts        <- wide2long((model)(times, par))
  datapts        <- do.call(rbind, lapply(1:length(times), function(i) datapts[intersect(which(datapts$time==times[i]),which(datapts$name==obs_names[i])),] ))
  return(list(sens, sens_inv, datapts))
}
###
trust_new <- function (parinit, rinit, rmax, parscale, iterlim = 100, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps), minimize = TRUE, blather = FALSE,
                       data, sigma, names, model, timeframe = c(0,50,.1), datainit, thresh = 1e-3, exact = FALSE, ...){
  
  datapts      <- wide2long((model)(datainit$time, parinit))
  datapts      <- do.call(rbind, lapply(1:length(datainit$time), function(i) datapts[intersect(which(datapts$time==datainit$time[i]),which(as.character(datapts$name)==as.character(datainit$name[i]))),] ))
  
  trafo <- Pdata(times = datapts$time, obs_names = datapts$name, par = parinit, model = model, thresh = thresh)
  
  parinit_data <- structure( as.numeric(datapts$value), names = colnames(trafo[[2]])) 
  attr(parinit,"deriv") <- trafo[[2]] #append derivs
  
  objfun <- normL2(data = data, model)
  
  if (!is.numeric(parinit)){
    stop("parinit not numeric")
  }
  if (!all(is.finite(parinit))){
    stop("parinit not all finite")
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
  if (!is.logical(minimize)){
    stop("minimize not logical")
  }
  
  r              <- rinit
  theta          <- parinit
  theta_data     <- parinit_data
  theta_data_act <- parinit_data
  out <- try(objfun(theta, ...))
  
  if (inherits(out, "try-error")) {
    warning("error in first call to objfun")
    return(list(error = out, argument = theta, converged = FALSE, 
                iterations = 0))
  }
  
  check.objfun.output(out, minimize, d) #checks for datatypes of input etc.
  
  if (!is.finite(out$value)){
    stop("parinit not feasible")
  }
  
  accept <- TRUE
  
  if (blather){
    theta.blather          <- NULL
    theta_data.blather     <- NULL
    theta_data_act.blather <- NULL
    theta.try.blather      <- NULL
    theta.try_data.blather <- NULL
    type.blather           <- NULL
    accept.blather         <- NULL
    r.blather              <- NULL
    stepnorm.blather       <- NULL
    rho.blather            <- NULL
    val.blather            <- NULL
    val.try.blather        <- NULL
    preddiff.blather       <- NULL
  }
  
  ###### MAIN LOOP ######
  for (iiter in 1:iterlim) {
    if (blather) {
      theta.blather      <- rbind(theta.blather, theta)
      theta_data.blather <- rbind(theta_data.blather, theta_data)
      theta_data_act.blather <- rbind(theta_data_act.blather, theta_data_act)
      r.blather <- c(r.blather, r)
      if (accept)
        val.blather <- c(val.blather, out$value)
      else val.blather <- c(val.blather, out.value.save)
    }
    
    if (accept) {
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
      gq <- as.numeric(t(eout$vectors) %*% g) #multiply eigenvectors of B with gradient
    }
    is.newton <- FALSE
    if (all(eout$values > 0)) {
      ptry <- as.numeric(-eout$vectors %*% (gq/eout$values))
      if (norm(ptry) <= r) 
        is.newton <- TRUE # if B is Hessian it will perform  a Newton step -> step is minimizer of quadratic function
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
          uout <- uniroot(fred, c(beta.dn, beta.up))
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
    
    #atm ptry is dataparameter step
    ptry_data <- ptry
    
    theta_dummy <- theta
    attr(theta_dummy, "deriv") <- NULL
    
    if(exact){
      datatry <- datapts
      datatry$value <- theta_data + ptry_data
      datatry$sigma <- 10^-5
      
      datatry <- as.datalist(datatry)
      
      objtry <- normL2(datatry, model)
      fittry <- trust(objfun = objtry, parinit = theta_dummy, rinit = 0.01, rmax = 1)
      fittry.converged <- fittry$converged
    }
    else{fittry.converged <- FALSE}
    if(fittry.converged){
      ptry <- fittry$argument-theta_dummy
      print("CONVERGED")
    }
    else{
      ptry      <- as.numeric(trafo[[2]] %*% ptry_data)
      if(!exact){
        print("option: not exact")
      }
      else{print("not converged")}
    }
    ###############################
    preddiff <- sum(ptry_data * (g + as.numeric(B %*% ptry_data)/2))
    
    if (rescale) {
      theta.try <- theta + ptry/parscale #STEP
      theta.try_data <- theta_data + ptry_data/parscale
    }
    else {
      theta.try <- theta + ptry #STEP
      theta.try_data <- theta_data + ptry_data
    }
    out <- try(objfun(theta.try, ...))
    if (inherits(out, "try-error")) 
      break
    check.objfun.output(out, minimize, d)
    ftry <- out$value
    if (!minimize) 
      ftry <- (-ftry)
    rho <- (ftry - f)/preddiff
    if (ftry < Inf) {
      is.terminate <- abs(ftry - f) < fterm || abs(preddiff) < mterm #termination tolerance check
    }
    else {
      is.terminate <- FALSE
      rho <- (-Inf)
    }
    if (is.terminate) {
      if (ftry < f) {
        accept <- TRUE
        theta      <- theta.try #value is smaller at the proposed step, so the new parameter is accepted
        theta_data <- theta.try_data
      }
    }
    else {
      if (rho < 1/4) {
        accept <- FALSE
        r <- r/4
      }
      else {
        accept <- TRUE
        theta      <- theta.try
        theta_data <- theta.try_data
        if (rho > 3/4 && (!is.newton)) 
          r <- min(2 * r, rmax)
      }
    }
    if (blather) {
      theta.try.blather <- rbind(theta.try.blather, theta.try)
      theta.try_data.blather <- rbind(theta.try_data.blather, theta.try_data)
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
    
    #new sensitivities
    #return(list(theta, theta_data))
    attr(theta, "deriv") <- NULL
    trafo <- Pdata(times = datapts$time, obs_names = datapts$name, par = theta, model = model, thresh = thresh)
    attr(theta, "deriv") <- trafo[[2]]
    theta_data_act <- theta_data
    theta_data <- trafo[[3]]$value
    print(theta_data_act)
    print(theta_data)
    print("-----")
  }
  if (inherits(out, "try-error")) {
    out <- list(error = out, argument = theta.try, converged = FALSE)
  }
  else {
    out <- try(objfun(theta, ...))
    if (inherits(out, "try-error")) {
      out <- list(error = out)
      warning("error in last call to objfun")
    }
    else {
      check.objfun.output(out, minimize, d)
    }
    out$argument <- theta
    out$converged <- is.terminate
  }
  out$iterations <- iiter
  if (blather) {
    dimnames(theta.blather)           <- NULL
    dimnames(theta_data.blather)      <- NULL
    dimnames(theta_data_act.blather)  <- NULL
    out$argpath                       <- theta.blather
    out$dataargpath                   <- theta_data.blather
    out$dataargpathact                <- theta_data_act.blather
    dimnames(theta.try.blather)       <- NULL
    dimnames(theta.try_data.blather)  <- NULL
    out$argtry                        <- theta.try.blather
    out$dataargtry                    <- theta.try_data.blather
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

#######################################################################################################################################################################

perfCheck <- function(parameters, timeframe, n_data, model, sigma, datainit, rinit, rmax, thresh, iterlim){
  do.call(rbind, mclapply(1:1000, function(i){
    off = structure(rnorm(length(parameters),0,1), names = parameters)
    parinit = structure(rnorm(length(parameters),0,1), names = parameters)
    par_data = parinit + off
    
    times = seq(from = timeframe[1], to = timeframe[2], by = timeframe[3])
    initmodel  <- wide2long((model)(times, parinit))
    
    #generate data
    data_time  <- sort(c(0, runif(n_data,min = 0, max = 50)))
    data       <- wide2long((model)(data_time, par_data))
    data$value <- rnorm(length(data$value), mean = data$value, sd = sigma[1])
    data$sigma <- sigma[1]
    data_list  <- as.datalist(data)
    
    #new trust
    start             <- proc.time()[3]
    dummy             <- try(trust_new(parinit = parinit, datainit = datainit, rinit = rinit, rmax = rmax, data = data_list, sigma = sigma, names = parameters, model = model, thresh = thresh, blather=TRUE, iterlim = iterlim))
    time_dat          <- proc.time()[3] - start
    if(!(inherits(dummy, "try-error")) ){fit_list          <- rbind(fit_list, data.frame(value = c(dummy$value, time_dat, dummy$iterations, dummy$converged), name = c("value", "time", "iterations", "converged"), optimizer = "dat"))}
    
    start                <- proc.time()[3]
    dummy                <- try(trust_new(parinit = parinit, datainit = datainit, rinit = rinit, rmax = rmax, data = data_list, sigma = sigma, names = parameters, model = g*x*p, thresh = thresh, blather=TRUE, iterlim = iterlim, exact = TRUE))
    time_dat_ex          <- proc.time()[3] - start
    if(!(inherits(dummy, "try-error")) ){fit_list             <- rbind(fit_list, data.frame(value = c(dummy$value, time_dat_ex, dummy$iterations, dummy$converged), name = c("value", "time", "iterations", "converged"), optimizer = "dat_ex"))}
    
    start                <- proc.time()[3]
    obj                  <- normL2(data_list, g*x*p)
    dummy                <- try(trust(obj, parinit = parinit, rinit = rinit, rmax = rmax, blather = TRUE, iterlim = iterlim))
    time_normal          <- proc.time()[3] - start
    if(!(inherits(dummy, "try-error")) ){fit_list             <- rbind(fit_list, data.frame(value = c(dummy$value, time_normal, dummy$iterations, dummy$converged), name = c("value", "time", "iterations", "converged"), optimizer = "normal"))}
    return(fit_list)
  }, mc.cores = 24, mc.preschedule = FALSE))
}