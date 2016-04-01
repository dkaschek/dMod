#' Return some useful forcing functions as strings
#' 
#' @param type Which function to be returned
#' @param parameters Named vector, character or numeric. Replace parameters by the corresponding valus
#' in \code{parameters}.
#' @return String with the function
#' @export
forcingsSymb <- function(type =c("Gauss", "Fermi", "1-Fermi", "MM", "Signal"), parameters = NULL) {
  
  type <- match.arg(type)
  fun <- switch(type,
                "Gauss"   = "(scale*exp(-(time-mu)^2/(2*tau^2))/(tau*2.506628))",
                "Fermi"   = "(scale/(exp((time-mu)/tau)+1))",
                "1-Fermi" = "(scale*exp((time-mu)/tau)/(exp((time-mu)/tau)+1))",
                "MM"      = "(slope*time/(1 + slope*time/vmax))",
                "Signal"  = "max1*max2*(1-exp(-time/tau1))*exp(-time*tau2)"
  )
  
  if(!is.null(parameters)) {
    fun <- replaceSymbols(names(parameters), parameters, fun)
  }
  
  return(fun)
  
}


#' Get coefficients from a character
#' 
#' @param char character, e.g. "2*x + y"
#' @param symbol single character, e.g. "x" or "y"
#' @return numeric vector with the coefficients
getCoefficients <- function(char, symbol) {
  
  pdata <- getParseData(parse(text = char, keep.source = TRUE))
  pdata <- pdata[pdata$terminal == TRUE, ] #  subset(pdata, terminal == TRUE)
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
#' @export
resolveRecurrence <- function (variables) {
  if(length(variables) > 1) {
    for (i in 1:(length(variables) - 1)) {
      newvariables <- c(variables[1:i], 
                        unlist(replaceSymbols(names(variables)[i],
                                              paste("(", variables[i], ")", sep = ""), 
                                              variables[(i + 1):length(variables)])))
      names(newvariables) <- names(variables)
      variables <- newvariables
    }
  }
  
  return(variables)
}

