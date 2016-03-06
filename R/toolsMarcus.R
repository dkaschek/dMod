

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
simulateData <- function(timesD, pouter, vars=observables, relE = 0.05, absE = 0.001){        

  ## Die Funktion setzt voraus, dass eine prediction function "x" im workspace ist :-(
  
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



#' Calculate analytical steady states
#' 
#' @param model Name of the file which the stoechiometric matric in csv format.
#' @param file Name of the file to which the steady-state equations are saved.
#'   Read this file with \code{\link{readRDS}}.
#' @param forcings Character vector with the names of the forcings
#' @param neglect Character vector with the names of species which must not be 
#'   used for solving the steady state equations.
#' @param outputFormat Define the output format. By default "R" generating dMod 
#'   compatible output. To obtain an output appropriate for d2d [1] "M" must be 
#'   selected.
#'   
#' @return Character vector of steady-state equations.
#'   
#' @references [1]
#' \url{https://bitbucket.org/d2d-development/d2d-software/wiki/Home}
#' 
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
steadyStates <- function(model, file, forcings = "", neglect = "", outputFormat = "R", sparsifyLevel = 2) {
  
  # Check if file is valid
  if (!is.character(file)) stop("File name must be specified")
  
  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {
    
    write.eqnlist(model, file = paste0(file, "_model.csv"))
    model <- paste0(file, "_model.csv")
    
  }
  
  # Calculate steady states.
  python.version.request("2.7")
  python.load(system.file("code/steadyStates.py", package = "dMod"))
  m_ss <- python.call("ODESS", model, forcings, neglect, outputFormat, sparsifyLevel)
  
  # Write steady states to disk.
  m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
    out <- eq[2]
    names(out) <- eq[1]
    return(out)
    }))
  saveRDS(object = m_ssChar, file = file)
}
