#' Simulate data for the currently loaded model
#'
#' @param timesD Desired time points for the data
#' @param pouter Outer parameter vector for which simulation is performed
#' @param vars Vector of strings of observables for which data shall be generated
#' @param relE Relative Error
#' @param absE Absolute Error
#'
#' @return List of data points for each condition
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#' 
#' @export
simulateData <- function(timesD, pouter, vars=observables, relE = 0.05, absE = 0.001){        
  pred <- x(timesD, pouter)
  out <- lapply(conditions, function(con){
    mydata <- wide2long(pred[[con]])
    sigma <- relE*mydata$value + absE
    mydata$value <- mydata$value + rnorm(dim(mydata)[1], 0, sigma)
    mydata <- cbind(mydata, sigma=sigma)
    subset(mydata, name%in%vars)
  }); names(out) <- conditions
  
  as.datalist(out)
  
}