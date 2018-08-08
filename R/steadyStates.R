#' Calculate analytical steady states
#' 
#' @param model Either name of the csv-file or the eqnlist of the model. If NULL, specify smatrix, states and rates by hand.
#' @param file Name of the file to which the steady-state equations are saved.
#'   Read this file with \code{\link{readRDS}}.
#' @param smatrix Numeric matrix, stiochiometry matrix of the system 
#' @param states Character vector, state vector of the system
#' @param rates Character vector, flux vector of the system
#' @param forcings Character vector with the names of the forcings
#' @param givenCQs Character vector with conserved quantities. Use the format c("A + pA = totA", "B + pB = totB"). If NULL, conserved quantities are automatically calculated.
#' @param neglect Character vector with names of states and parameters that must not be used for solving the steady-state equations
#' @param sparsifyLevel numeric, Upper bound for length of linear combinations used for simplifying the stoichiometric matrix
#' @param outputFormat Define the output format. By default "R" generating dMod 
#'   compatible output. To obtain an output appropriate for d2d [1] "M" must be 
#'   selected.
#'   
#' @return Character vector of steady-state equations.
#'   
#' @references [1]
#' \url{https://bitbucket.org/d2d-development/d2d-software/wiki/Home}
#' 
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'   
#' @export
#' @importFrom utils write.table
steadyStates <- function(model, file=NULL, smatrix = NULL, states = NULL, rates = NULL, forcings = NULL, givenCQs = NULL, neglect=NULL, sparsifyLevel = 2, outputFormat = "R") {
  
  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {    
    write.eqnlist(model, file = paste0(file, "_model.csv"))
    model <- paste0(file, "_model.csv")    
  }
  if(!is.null(smatrix)){
    write.table(smatrix, file="smatrix.csv", sep = ",")
    smatrix=TRUE
  }
  
  # Calculate steady states.
  python.version.request("2.7")  
  rPython::python.load(system.file("code/steadyStates.py", package = "dMod"))
  m_ss <- rPython::python.call("ODESS", model, smatrix, as.list(states), as.list(rates), as.list(forcings), as.list(givenCQs), as.list(neglect), sparsifyLevel, outputFormat)
  
  # Write steady states to disk.
  if(length(m_ss)>1){    
    m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
      out <- eq[2]
      names(out) <- eq[1]
      return(out)
    }))
    if(!is.null(file) & is.character(file)){
      saveRDS(object = m_ssChar, file = file)
    }
  } else return(0)
}

