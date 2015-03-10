#' Compare data in integration output to compute residuals
#' 
#' @param data data.frame with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param out output of ode(), optionally augmented with attributes 
#' "deriv" (output of ode() for the sensitivity equations) and
#' "parameters" (character vector of parameter names, a subsest of those 
#' contained in the sensitivity equations). If "deriv" is given, also "parameters"
#' needs to be given.
#' @return data.frame with the original data augmented by columns "prediction" (
#' numeric, the model prediction), "residual" (numeric, difference between
#' data value and prediction), "weighted.residual" (numeric, residual devided
#' by sigma). If "deriv" was given, the returned data.frame has an 
#' attribute "deriv" (data.frame with the derivatives of the residuals with 
#' respect to the parameters).
res <- function (data, out) {
  
  times <- sort(unique(data$time))
  names <- as.character(unique(data$name))
  match.times <- match(data$time, times)
  match.names <- match(data$name, names)
  match.coords <- cbind(match.times, match.names)
  index <- apply(match.coords, 1, function(m) length(times)*(m[2]-1) + m[1])
  
  
  
  names.extended <- rep(names, each=length(times))
  values.extended <- matrix(NA, nrow = length(times)*length(names), ncol=3, 
                            dimnames = list(NULL, c("time", "value", "sigma")))
  values.extended[index, "time"] <- data$time
  values.extended[index, "value"] <- data$value
  values.extended[index, "sigma"] <- data$sigma
  
  data <- data.frame(name = names.extended, values.extended)
  
  
  outtimes <- out[, 1]
  index <- match(times, outtimes)
  nout <- as.vector(out[index, names])
  deriv <- attr(out, "deriv")
  pars <- attr(out, "parameters")
  sensData <- NULL
  if (!is.null(deriv)) {
    sensnames <-as.vector(outer(names, pars, paste, sep="."))
    noutSens <- as.vector(deriv[index, sensnames])
    M <- matrix(noutSens, ncol = length(pars), dimnames = list(NULL, pars))
    sensData <- data.frame(time = data$time, name = data$name, M)
  }
  residuals <- data$value - nout
  residuals[is.na(residuals)] <- 0
  weighted.residuals <- (data$value - nout)/data$sigma
  weighted.residuals[is.na(weighted.residuals)] <- 0
  data <- cbind(data, prediction = nout, residual = residuals, 
                weighted.residual = weighted.residuals)
  data <- data[c("time", "name", "value", "prediction", "sigma", 
                 "residual", "weighted.residual")]
  attr(data, "deriv") <- sensData
  return(data)
}

#' Compute the weighted residual sum of squares
#' 
#' @param nout data.frame (result of \link{res})
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
wrss <- function(nout) {
  
  obj <- sum(nout$weighted.residual^2)
  grad <- NULL
  hessian <- NULL
  
  if(!is.null(attr(nout, "deriv"))) {
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
  
    sens <- as.matrix(attr(nout, "deriv")[,-(1:2)])
    grad <- -as.vector(2*matrix(nout$residual/nout$sigma^2, nrow=1)%*%sens)
    names(grad) <- colnames(sens)
    hessian <- 2*t(sens/nout$sigma)%*%(sens/nout$sigma)
    
    
  }
  
  out <- list(value=obj, gradient=grad, hessian=hessian)
  class(out) <- c("obj", "list")
  
  return(out)

}






#' Translate data into time-continuous data representations
#' 
#' @param data \code{data.frame} with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param tau Numeric. The uncertainty of a single data point is smeared by a Gaussian with variance tau².
#' @param ngrid The number of time points returned.
#' @param type Character, the type of interpolation, i.e. "linear", "spline", "logspline" or "smooth".
#' @return a \code{data.frame} with columns "name", "time" and "value".
data2forc <- function(data, tau=NULL, ngrid = 1e3, type="linear") {
  
  obs <- unique(data$name)
  obsD <- paste0(obs, "D")
  weightsD <- paste0("weight", obs, "D")
  tmin <- min(data$time)
  tmax <- max(data$time)
  times <- seq(tmin, tmax, len=ngrid)
  
  # data points
  dataout <- do.call(rbind, lapply(obs, function(o) {
    
    subdata <- subset(data, name == o)
    if(type=="linear") {
      t <- subdata$time
      x <- subdata$value
      xoft <- approxfun(t, x)
      xvalues <- xoft(times)
      
      xofts <- splinefun(times, xvalues)
      
      out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xvalues),
                   data.frame(name = paste0(o, "DDot"),time = times, value = xofts(times, deriv=1)))
    }
    if(type=="spline") {
      t <- subdata$time
      x <- subdata$value
      xoft <- splinefun(t, x)
      out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xoft(times)),
                   data.frame(name = paste0(o, "DDot"), time = times, value = xoft(times, deriv=1)))
    }
    if(type=="logspline") {
      t <- subdata$time
      
      if(any(subdata$value <= 0)) {
        # no log
        x <- subdata$value 
      } else {
        # log
        subdata$value[subdata$value < 1e-5*max(subdata$value)] <- 1e-5*max(subdata$value)
        x <- log(subdata$value)
      }
      
      xoft <- splinefun(t, x)
      
      if(any(subdata$value <= 0)) {
        # no log
        out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = xoft(times)),
                     data.frame(name = paste0(o, "DDot"), time = times, value = splinefun(times, xoft(times))(times, deriv=1)))
      } else {
        # log
        out <- rbind(data.frame(name = paste0(o, "D"),time = times, value = exp(xoft(times))),
                     data.frame(name = paste0(o, "DDot"), time = times, value = splinefun(times, exp(xoft(times)))(times, deriv=1)))
      }
      
      
    }
    if (type == "smooth") {
      t <- subdata$time
      x <- subdata$value
      xoft <- function(newt, ...) predict(smooth.spline(t, x, df = 5), newt, ...)$y
      out <- rbind(data.frame(name = paste0(o, "D"), time = times, 
                              value = xoft(times)), data.frame(name = paste0(o, 
                                                                             "DDot"), time = times, value = xoft(times, deriv = 1)))
      
    }
    
    return(out)
    
  }))
  
  # sigma value  
  
  
  
  if(is.null(tau)) {
    Deltat <- (tmax-tmin)/(dim(data)[1]/length(obs))
    tau <- Deltat/5
  }
  
  if(tau == 0) {
    weightout <- do.call(rbind, lapply(obs, function(o) {
      
      subdata <- subset(data, name == o & !is.na(value) & !is.na(sigma))
      
      t <- subdata$time
      x <- log(1/subdata$sigma^2)
      xoft <- splinefun(t, x)
      
      
      rbind(data.frame(name = paste0("weight", o, "D"), time = times, value = exp(xoft(times))),
            data.frame(name = paste0("weight", o, "DDot"), time=times, value = splinefun(times, exp(xoft(times)))(times, deriv=1)))
      
    }))
    
  } else {
    weightout <- do.call(rbind, lapply(obs, function(o) {
      
      subdata <- subset(data, name == o & !is.na(value) & !is.na(sigma))
      support <- subdata$time
      supfactor <- rep(1, length(support))
      if(any(support==tmin)) supfactor[support==tmin] <- 2
      if(any(support==tmax)) supfactor[support==tmax] <- 2
      sigma <- subdata$sigma
      weights <- Reduce("+", lapply(1:length(support), function(i) supfactor[i]*dnorm(times, support[i], tau) * (1/sigma[i]^2)))
      weightsdot <- splinefun(times, weights)(times, deriv=1)
      rbind(data.frame(name = paste0("weight", o, "D"), time = times, value = weights),
            data.frame(name = paste0("weight", o, "DDot"), time = times, value = weightsdot))
    }))
    
  }
  
  out <- rbind(dataout, weightout)
  
  return(out)
  
  
}