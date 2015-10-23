



## Methods for the class parvec ------------------------------------------------


#' Pretty printing of parameter transformations, class par
#' 
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
print.parvec <- function(p, ...) {
  
  ## Diese Funktion mergen mit print.eqnvec
  ## Für parvec neue Funktion schreiben
  
  
#   # Assemble parameters
#   hInner <- "Inner"
#   hOuter <- "Outer"
#   arglist <- list(...)
#   if("linewidth" %in% names(arglist)) linewidth <- arglist$linewidth else linewidth <- 79
#   
#   maxNameWidth <- max(nchar(names(equations)), nchar(hInner))  
#   equationWidth <- if (linewidth - maxNameWidth - 3 > 9) linewidth - maxNameWidth -3
#   else 10
#   
#   
#   # Assemble and print parameter transformation table
#   cat("Table of parameter transformations\n")
#   for (i in seq(1, length(equations))) {
#     eq <- strelide(equations[i], equationWidth, where = "right")
#     eqName <- strpad(names(equations[i]), maxNameWidth, where = "left")
#     cat(eqName, " = ", eq, "\n", sep = "")
#     if (!(i %% 10)) {
#       cat("\n")
#     }
#   }
}


#' @export
"[.parvec" <- function(x, ...) {
  out <- unclass(x)[...]
  deriv <- attr(x, "deriv")[..., ]
  parvec(out, deriv = deriv)
  return(out)
}



