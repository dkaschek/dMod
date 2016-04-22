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

