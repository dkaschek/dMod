

#' Check which Python versions are installed on the system
#' 
#' @param version NULL or character. Check for specific version
#' @return Character vector with the python versions and where they are located.
#' @export
python.version.sys <- function(version = NULL) {
  
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
python.version.rpython <- function() {
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
python.version.request <- function(version) {
  
  # Is rPythen installed and linked against requested python version?
  m_installed <- "rPython" %in% installed.packages()[, 1]
  if (m_installed) {
    m_curVersion <- python.version.rpython()
    if (m_curVersion == version)
      return(TRUE)
  }
  
  # rPython not installed or linked agains wrong python version.
  m_sysVersion <- python.version.sys(version)
  if (is.null(m_sysVersion)) {
    msg <- paste0("Requested python version ", version, " not available\n",
                  "Your options are\n")
    cat(msg)
    print(python.version.sys())
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
    install.packages("rPython")
  }
}
