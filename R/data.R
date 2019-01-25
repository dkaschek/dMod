#' Compare data and model prediction by computing residuals
#' 
#' @param data data.frame with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param out output of ode(), optionally augmented with attributes 
#' "deriv" (output of ode() for the sensitivity equations) and
#' "parameters" (character vector of parameter names, a subsest of those 
#' contained in the sensitivity equations). If "deriv" is given, also "parameters"
#' needs to be given.
#' @param err output of the error model function
#' @return data.frame with the original data augmented by columns "prediction" (
#' numeric, the model prediction), "residual" (numeric, difference between
#' prediction and data value), "weighted.residual" (numeric, residual devided
#' by sigma). If "deriv" was given, the returned data.frame has an 
#' attribute "deriv" (data.frame with the derivatives of the residuals with 
#' respect to the parameters).
#' @export
#' @import cOde
#' @importFrom stats setNames
res <- function(data, out, err = NULL) {
  
  data$name <- as.character(data$name)
  
  # Unique times, names and parameter names
  times <- sort(unique(data$time))
  names <- unique(data$name)
  
  # Match data times/names in unique times/names
  data.time <- match.num(data$time, times)
  data.name <- match(data$name, names)
  
  # Match unique times/names in out times/names
  time.out <- match.num(times, out[,1])
  name.out <- match(names, colnames(out))
  

  # Match data times/names in out times/names
  timeIndex <- time.out[data.time]
  nameIndex <- name.out[data.name]
  prediction <- sapply(1:nrow(data), function(i) out[timeIndex[i], nameIndex[i]]) 

  # Propagate derivatives if available
  deriv <- attr(out, "deriv")
  deriv.data <- NULL    
    
  # Propagate derivatives of err model if available
  deriv.err <- attr(err, "deriv")
  deriv.err.data <- NULL
  
  # Set value to loq if below loq
  data$value <- pmax(data$value, data$lloq)
  is.bloq <- data$value <= data$lloq
  
  
  
  if (!is.null(deriv)) {
  
    pars <- unique(unlist(lapply(strsplit(colnames(deriv)[-1], split = ".", fixed = TRUE), function(i) i[2])))
    sensnames <- as.vector(outer(names, pars, paste, sep = "."))
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
  
 
  # Modifications if error model is available
  if (!is.null(err)) {
    
    time.err <- match.num(times, err[,1])
    name.err <- match(names, colnames(err))
    timeIndex <- time.err[data.time]
    nameIndex <- name.err[data.name]
    errprediction <- sapply(1:nrow(data), function(i) err[timeIndex[i], nameIndex[i]]) 
    data$sigma[!is.na(errprediction)] <- errprediction[!is.na(errprediction)]

    
    if (!is.null(deriv.err)) {
      
      pars <- unique(unlist(lapply(strsplit(colnames(deriv.err)[-1], split = ".", fixed = TRUE), function(i) i[2])))
      sensnames <- as.vector(outer(names, pars, paste, sep = "."))
      # Match names to the corresponding sensitivities in sensnames
      names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)))
      # Get positions of sensnames in colnames of deriv
      sensnames.deriv <- match(sensnames, colnames(deriv.err))
      # Get the columns in deriv corresponding to data names
      derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
      # Derivatives of the prediction
      deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) submatrix(deriv.err, timeIndex[i], derivnameIndex[, i])))
      colnames(deriv.prediction) <- pars
      deriv.prediction[is.na(deriv.prediction)] <- 0
      
      deriv.err.data <- data.frame(time = data$time, name = data$name, deriv.prediction)
      
    }
    
    
  }
  
  
  # Compute residuals
  residuals <- prediction - data$value 
  weighted.residuals <- (prediction - data$value)/data$sigma
  
  data[["prediction"]] <- prediction
  data[["residual"]] <- residuals
  data[["weighted.residual"]] <- weighted.residuals
  data[["bloq"]] <- is.bloq
  
  objframe(data, deriv = deriv.data, deriv.err = deriv.err.data)
  
}


#' Time-course data for the JAK-STAT cell signaling pathway
#'
#' Phosphorylated Epo receptor (pEpoR), phosphorylated STAT in the
#' cytoplasm (tpSTAT) and total STAT (tSTAT) in the cytoplasmhave been 
#' measured at times 0, ..., 60.
#'
#' @name jakstat
#' @docType data
#' @keywords data
NULL


# Match with numeric tolerance 
match.num <- function(x, y, tol = 1e-8) {
  
  digits <- -log10(tol)
  match(round(x, digits), round(y, digits))
  
} 

