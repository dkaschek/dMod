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
    names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)))
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv))
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
    # Derivatives of the prediction
    deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) submatrix(deriv, timeIndex[i], derivnameIndex[, i])))
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
  #attr(data, "deriv") <- deriv.data
  
  objframe(data, deriv = deriv.data)
  
}






