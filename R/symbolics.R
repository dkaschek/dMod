#' Get coefficients from a character
#' 
#' @param char character, e.g. "2*x + y"
#' @param symbol single character, e.g. "x" or "y"
#' @return numeric vector with the coefficients
#' @examples getCoefficients("2*x + x + y", "x")
getCoefficients <- function(char, symbol) {
  
  pdata <- getParseData(parse(text=char))
  pdata <- subset(pdata, terminal==TRUE)
  symbolPos <- which(pdata$text == symbol)
  coefficients <- rep(1, length(symbolPos))
  
  hasCoefficient <- rep(FALSE, length(symbolPos))
  hasCoefficient[symbolPos > 1] <- (pdata$text[symbolPos[symbolPos > 1] - 1] == "*")
  coefficients[hasCoefficient] <- pdata$text[symbolPos[hasCoefficient]-2]
  
  return(as.numeric(coefficients))
  
  
  
  
  
}


#' Place top elements into bottom elemens
#' 
#' @param variables named character vector
#' @details If the names of top vector elements occur in the bottom of the vector, 
#' they are replaced by the character of the top entry. Useful for steady state conditions.
#' @return named character vector of the same length as \code{variables}
#' @examples resolveRecurrence(c(A = "k1*B/k2", C = "A*k3+k4", D="A*C*k5"))
resolveRecurrence <- function (variables) 
{
  for (i in 1:(length(variables) - 1)) {
    newvariables <- c(variables[1:i], 
                      unlist(replaceSymbols(names(variables)[i],
                                            paste("(", variables[i], ")", sep = ""), 
                                            variables[(i + 1):length(variables)])))
    names(newvariables) <- names(variables)
    variables <- newvariables
  }
  return(variables)
}




# sensitivitiesSymb.eqn <- function(f, inputs = NULL) {
#   
#   v <- attr(f, "rates"); names(v) <- paste0("v", 1:length(v))
#   S <- attr(f, "SMatrix")
#   S[is.na(S)] <- 0
#   states <- attr(f, "species")
#   parameters <- getSymbols(v, exclude=c(states, inputs))
#   
#   
#   
#   Dxv <- jacobianSymb(v, states)
#   Dpv <- jacobianSymb(v, parameters)
#   Dxv.zero <- names(Dxv)[Dxv == "0"]
#   Dpv.zero <- names(Dpv)[Dpv == "0"]
#   
#   x.x <- matrix(apply(expand.grid.alt(states, states), 1, paste, collapse = "."), length(states))
#   x.p <- matrix(apply(expand.grid.alt(states, parameters), 1, paste, collapse = "."), length(states))
#   v.x <- matrix(apply(expand.grid.alt(names(v), states), 1, paste, collapse = "."), length(v))
#   v.p <- matrix(apply(expand.grid.alt(names(v), parameters), 1, paste, collapse = "."), length(v))
#   
#   v.x[match(Dxv.zero, v.x)] <- "0"
#   v.p[match(Dpv.zero, v.p)] <- "0"
#   
#   
#   #print(prodSymb(v.x, x.x))
#   #print(S)
#   x.x.t <- structure(prodSymb(t(S), prodSymb(v.x, x.x)), names = x.x)
#   
#   x.p.t <- structure(prodSymb(t(S), sumSymb(prodSymb(v.x, x.p), v.p)), names = x.p)
#   
#   
#   rates <- c(Dxv[!names(Dxv)%in%Dxv.zero], 
#              Dpv[!names(Dpv)%in%Dpv.zero])
#   
#   out <- c(x.x.t, 
#            x.p.t)
#   
#   attr(out, "rates") <- rates
#   
#   
#   return(out)
#   
#   
#   
# }
