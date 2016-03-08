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
#' @param sparsifyLevel numeric, set the level to which the stoichiometric matrix is
#'   simplified by searching for clever linear combinations.   
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
  rPython::python.load(system.file("code/steadyStates.py", package = "dMod"))
  m_ss <- rPython::python.call("ODESS", model, forcings, neglect, outputFormat, sparsifyLevel)
  
  # Write steady states to disk.
  m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
    out <- eq[2]
    names(out) <- eq[1]
    return(out)
    }))
  saveRDS(object = m_ssChar, file = file)
}
