#' Compare data and model prediction by computing residuals
#' 
#' @param data data.frame with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param out output of ode(), optionally augmented with attributes 
#' "deriv" (output of ode() for the sensitivity equations) and
#' "parameters" (character vector of parameter names, a subsest of those 
#' contained in the sensitivity equations). If "deriv" is given, also "parameters"
#' needs to be given.
#' @return data.frame with the original data augmented by columns "prediction" (
#' numeric, the model prediction), "residual" (numeric, difference between
#' prediction and data value), "weighted.residual" (numeric, residual devided
#' by sigma). If "deriv" was given, the returned data.frame has an 
#' attribute "deriv" (data.frame with the derivatives of the residuals with 
#' respect to the parameters).
#' @export
#' @import cOde
res <- function (data, out) {
  
  # Unique times, names and parameter names
  times <- sort(unique(data$time))
  names <- as.character(unique(data$name))
  pars <- attr(out, "parameters")
  
  # Match data times/names in unique times/names
  data.time <- match(data$time, times)
  data.name <- match(data$name, names)
  
  # Match unique times/names in out times/names
  time.out <- match(times, out[,1])
  name.out <- match(names, colnames(out))
  
  # Match data times/names in out times/names
  timeIndex <- time.out[data.time]
  nameIndex <- name.out[data.name]
  prediction <- sapply(1:nrow(data), function(i) out[timeIndex[i], nameIndex[i]]) 
  
  # Propagate derivatives if available
  deriv <- attr(out, "deriv")
  deriv.data <- NULL
  if (!is.null(deriv)) {
    sensnames <- as.vector(outer(names, pars, paste, sep="."))
    # Match names to the corresponding sensitivities in sensnames
    names.sensnames <- apply(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)), 1, identity)
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv))
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
    # Derivatives of the prediction
    deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) deriv[timeIndex[i], derivnameIndex[, i]]))
    colnames(deriv.prediction) <- pars
    
    deriv.data <- data.frame(time = data$time, name = data$name, deriv.prediction)
  }
  
  # Compute residuals
  residuals <- prediction - data$value 
  weighted.residuals <- (prediction - data$value)/data$sigma
  data <- cbind(data, prediction = prediction, residual = residuals, 
                weighted.residual = weighted.residuals)
  data <- data[c("time", "name", "value", "prediction", "sigma", 
                 "residual", "weighted.residual")]
  attr(data, "deriv") <- deriv.data
  return(data)
}


#' Compute the weighted residual sum of squares
#' 
#' @param nout data.frame (result of \link{res})
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
  
  out <- list(value=obj, gradient=grad, hessian=hessian)
  class(out) <- c("obj", "list")
  
  return(out)

}


#' Generate dummy list of class \code{obj} from named numeric
#' 
#' @param p Names numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.obj(p)
#' @export
as.obj <- function(p) {
  
  obj <- list(
    value = 0,
    gradient = structure(rep(0, length(p)), names = names(p)),
    hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
  class(obj) <- "obj"
  
  return(obj)
  
}





#' Translate data into time-continuous data representations
#' 
#' @param data \code{data.frame} with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param tau Numeric. The uncertainty of a single data point is smeared by a Gaussian with variance tau².
#' @param ngrid The number of time points returned.
#' @param type Character, the type of interpolation, i.e. "linear", "spline", "logspline" or "smooth".
#' @return a \code{data.frame} with columns "name", "time" and "value".
#' @export
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