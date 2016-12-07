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
  
  myframe <- parframe(
    cbind(
      index = 1:nrow(sample.matrix),
      as.data.frame(sample.matrix)
    ),
    metanames = "index",
    parameters= colnames(sample.matrix)
  )
  
  return(myframe)

}

#' Reparameterization
#' 
#' @param expr character of the form \code{"lhs ~ rhs"} where \code{rhs}
#' reparameterizes \code{lhs}. Both \code{lhs} and \code{rhs}
#' can contain a number of symbols whose values need to be passed by the \code{...} argument.
#' @param trafo character or equation vector where the replacement takes place
#' @param ... pass symbols as named arguments
#' @return an equation vector with the reparameterization.
#' @details Left and right-hand side of \code{expr} are searched for symbols. If separated by
#' "_", symbols are recognized as such, e.g. in \code{Delta_x} where the symbols are 
#' "Delta" and "x". Each symbol for which values (character or numbers) are passed by the
#' \code{...} argument is replaced.
#' @export
#' @examples
#' innerpars <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' 
#' trafo <- repar("x ~ x", x = innerpars)
#' trafo <- repar("x ~ y", trafo, x = names(constraints), y= constraints)
#' trafo <- repar("x ~ exp(x)", trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition", trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ...) {
 
  if (inherits(expr, "formula")) expr <- deparse(expr)
   
  parsed.expr <- as.character(as.formula(gsub("_", ":", expr, fixed = TRUE)))
  lhs <- parsed.expr[2]
  lhs.symbols <- getSymbols(lhs)
  rhs <- parsed.expr[3]
  rhs.symbols <- getSymbols(rhs)
  
  replacements <- as.data.frame(list(...), stringsAsFactors = FALSE)
  
  lhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], lhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  rhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], rhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  if (is.null(trafo)) trafo <- structure(lhs, names = lhs)
  trafo <- replaceSymbols(lhs, rhs, trafo)
  
  return(trafo)
  
  
}

