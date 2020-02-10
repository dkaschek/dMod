#' Calculate analytical steady states. 
#' 
#' @description This function follows the method published in [1]. The determined steady-state solution is tailored to parameter estimation. Please note that kinetic parameters might be fixed for solution of steady-state equations. Note that additional parameters might be introduced to ensure positivity of the solution.
#' @description The function calls a python script via rPython. Usage problems might occur when different python versions are used. The script was written and tested for python 2.7.12, sympy 0.7.6 and numpy 1.8.2.
#' @description Recently, users went into problems with RJSONIO when rPython was used. Unless a sound solution is available, please try to reinstall RJSONIO in these cases.
#' 
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
#'   compatible output. To obtain an output appropriate for d2d [2] "M" must be 
#'   selected.
#'   
#' @return Character vector of steady-state equations.
#'   
#' @references [1]
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863410/}
#' @references [2]
#' \url{https://github.com/Data2Dynamics/d2d}
#' 
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'   
#' @export
#' @importFrom utils write.table
#' @example inst/examples/steadystates.R
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
  python_version_request("2.7")  
  rPython::python.load(system.file("code/steadyStates.py", package = "dMod"))
  m_ss <- rPython::python.call("ODESS", model, smatrix, as.list(states), as.list(rates), as.list(forcings), as.list(givenCQs), as.list(neglect), sparsifyLevel, outputFormat)
  
  # Write steady states to disk.
  if(length(m_ss)>1){    
    m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
      out <- eq[2]
      names(out) <- eq[1]
      return(out)
    }))
    if(!is.null(file) & is.character(file))
      saveRDS(object = m_ssChar, file = file)
    
    return(m_ssChar)
  } else return(0)
}



#' Check which Python versions are installed on the system
#' 
#' @param version NULL or character. Check for specific version
#' @return Character vector with the python versions and where they are located.
#' @export
python_version_sys <- function(version = NULL) {
  
  # Which python versions are installed on the system
  m_sysPath <- strsplit(Sys.getenv("PATH"), ":")
  m_sysPath <- m_sysPath[[1]]
  m_python <- do.call(rbind, lapply(m_sysPath, function(p) {
    m_py <- dir(p, pattern = "^python[0-9,.]+$")
    if (length(m_py) != 0) {
      return(file.path(p, m_py))
    } else {
      return(NULL)
    }
  }))
  
  m_version <- strsplit(m_python, "python")
  m_version <- lapply(m_version, function(p) {
    return(p[2])
  })
  
  m_python <- as.data.frame(m_python)
  names(m_python) <- m_version
  
  
  if (is.null(version)) {
    return(m_python)  
  } else {
    # Is requested version available
    if (any(m_version == version)) {
      return(as.character(m_python[[version]]))
      attr(out, "version") <- m_version
    } else {
      return(NULL)
    }
  }
}


#' Get the Python version to which rPython is linked
#' 
#' @return The Python version and additional information
#' @export
python_version_rpython <- function() {
  rPython::python.exec(c("def ver():", "\timport sys; return list(sys.version_info)"))
  m_info <- as.data.frame(rPython::python.call("ver"))
  names(m_info) <- c("major", "minor", "micro", "releselevel", "serial")
  
  m_version <- paste0(m_info[[1]], ".", m_info[[2]])
  attr(m_version, "info") <- m_info
  
  return(m_version)
}


#' Check if rPython comes with the correct Python version
#' 
#' @description rPython is liked against a certain Python version found on the system.
#' If Python code called from R requires a specific Python version, the rPython package
#' needs to be reinstalled. This functions helps to do this in one line.
#' 
#' @param version character indicating the requested Python version
#' 
#' @return TRUE if rPython is linked against the requested version. Otherwise, the user
#' is asked if rPython should be reinstalled with the correctly linked Python version.
#' @export
python_version_request <- function(version) {
  
  # Is rPythen installed and linked against requested python version?
  m_installed <- "rPython" %in% utils::installed.packages()[, 1]
  if (m_installed) {
    m_curVersion <- python_version_rpython()
    if (m_curVersion == version)
      return(TRUE)
  }
  
  # rPython not installed or linked agains wrong python version.
  m_sysVersion <- python_version_sys(version)
  if (is.null(m_sysVersion)) {
    msg <- paste0("Requested python version ", version, " not available\n",
                  "Your options are\n")
    cat(msg)
    print(python_version_sys())
  } else {
    # If rPython is already installed, double check with user
    if (m_installed) {
      msg <- paste0("rPython is installed on your system using python ", m_curVersion, "\n",
                    "Proceeding the installation enables ", version, " for R\n",
                    "This will prevent python programms needing version ", m_curVersion, " from executing\n",
                    "Do you want to abort [1] or procede [2] with the installation?\n")
      cat(msg)
      m_go <- scan(nmax = 1, what = integer(), quiet = TRUE)
      if (m_go != 2) {
        stop("Installation aborted")
      }
    }
    
    try(detach("package:rPython", unload = TRUE), silent = TRUE)
    Sys.setenv(RPYTHON_PYTHON_VERSION = version)
    utils::install.packages("rPython")
  }
}

