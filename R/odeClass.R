## Methods of class odemodel

#' @export
print.odemodel <- function(x, ...) {
  
  func <- x$func
  extended <- x$extended
  
  suppressWarnings({
  
  cat("DLL name ODE: ", func, "\n", sep = "")
  cat("DLL name SENS: ", extended, "\n", sep = "")
  cat("Equations:\n", sep = "")
  print(as.eqnvec(attr(func, "equations")))
  cat("States:\n", sep = "")
  print(sort(attr(func, "variables")))
  cat("Parameters:\n", sep = "")
  print(sort(attr(func, "parameters")))
  cat("Forcings:\n", sep = "")
  print(sort(attr(func, "forcings")))
  
  })
}