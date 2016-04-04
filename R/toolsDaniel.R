#' Produce data frames for plotting
#' 
#' @param prediction object of class prediction list
#' @param data object of class data list
#' @param condition.grid data.frame with row.names according to the condition names.
#' @param substitute Substitute "_" by character \code{substitute}.
#' @param ... arguments going to subsetting of prediction and data
#' @return list with prediction (data.frame) and data (data.frame)
#' @export
ggdata <- function(prediction, data = NULL, condition.grid = attr(data, "condition.grid"), substitute = "_", ...) {
  
  prediction <- wide2long(prediction)
  data <- lbind(data)
  
  for (C in colnames(condition.grid)) {
    prediction[, C] <- condition.grid[as.character(prediction$condition), C]
    data[, C] <- condition.grid[as.character(data$condition), C]
  }
  
  prediction <- subset(prediction, ...)
  if (!is.null(data)) data <- subset(data, ...)
  
  n1 <- nrow(prediction)
  
  out <- combine(prediction, data)
  for (n in colnames(out)) {
    mylevels <- levels(out[, n]) 
    if (!is.null(mylevels)) levels(out[, n]) <- gsub("_", substitute, mylevels, fixed = TRUE)
  }
  
  return(list(prediction = out[1:n1,], data = out[-(1:n1),]))
  
}
#' Generate sample for multi-start fit
#' 
#' @param center named numeric, the center around we sample
#' @param samplefun character, indicating the random number generator,
#' defaults to \code{"rnorm"}.
#' @param fits length of the sample
#' @param ... arguments going to \code{samplefun}
#' @return matrix with the parameter samples
#' @export
mssample <- function(center, samplefun = "rnorm", fits = 20, ...) {
  
  sample.matrix <- do.call(rbind, lapply(1:fits, function(i) 
    do.call(samplefun, c(list(n = length(center), mean = center), list(...)))))
  
  colnames(sample.matrix) <- names(center)
  
  return(sample.matrix)

}

