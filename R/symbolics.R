#' Return some useful forcing functions as strings
#' 
#' @param type Which function to be returned
#' @param parameters Named vector, character or numeric. Replace parameters by the corresponding valus
#' in \code{parameters}.
#' @return String with the function
#' @export
forcingsSymb <- function(type =c("Gauss", "Fermi", "1-Fermi", "MM", "Signal", "Dose"), parameters = NULL) {
  
  type <- match.arg(type)
  
  # INPUT1 (differentiable box)
  #fn <- "(1/(1+exp(k*(time-T1))))*(exp(k*(time-T2))/(1+exp(k*(time-T2))))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  #fn <- "((exp(-k*(time-T1))/(exp(-k*(time-T1))+1))/(exp(-k*(time-T2))+1))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  fn1 <- "(.5*(1.-tanh(.5*k*(time-T1))))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  fn2 <- "(.5*(1.+tanh(.5*k*(time-T2))))" # T1 = start, T2 = end, k/4 = +-steepness in T1 and T2
  integral <- "(Tduration/(exp(20)-1))" # 100 = k*Tduration
  
  k <- "(20/Tduration)"
  T1 <- "(Tlag+Tinit)"
  T2 <- "(Tlag+Tinit+Tduration)"
  
  INPUT <- paste0("Dose*(", fn1, ")*(", fn2, ")/", integral)
  INPUT <- replaceSymbols(c("k", "T1", "T2"), c(k, T1, T2), INPUT)
  INPUT <- paste0("(", INPUT, ")")
  
  
  
  fun <- switch(type,
                "Gauss"   = "(scale*exp(-(time-mu)^2/(2*tau^2))/(tau*2.506628))",
                "Fermi"   = "(scale/(exp((time-mu)/tau)+1))",
                "1-Fermi" = "(scale*exp((time-mu)/tau)/(exp((time-mu)/tau)+1))",
                "MM"      = "(slope*time/(1 + slope*time/vmax))",
                "Signal"  = "(max1*max2*(1-exp(-time/tau1))*exp(-time*tau2))",
                "Dose"   = INPUT
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





#' Find integer-null space of matrix A 
#  this function is written along the lines of matlab's function null(A,'r'), where 'r' specifies that integer solutions are returned instead of orthonormal vectors
#' 
#' @param A matrix for which the null space is searched
#' @param tol tolerance to find pivots in rref-function below
#' @return null space of A with only integers in it
#' 
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
nullZ <- function(A, tol=sqrt(.Machine$double.eps)) {
  
  ret <- rref(A) # compute reduced row echelon form of A
  ret[[1]] -> R # matrix A in rref 
  ret[[2]] -> pivcol #columns in which a pivot was found
  
  n <- ncol(A) # number of columns of A
  r <- length(pivcol) # rank of reduced row echelon form
  nopiv <- 1:n
  nopiv <- nopiv[-pivcol]  # columns in which no pivot was found
  
  Z <- mat <- matrix(0, nrow = n, ncol = n-r) # matrix containing the vectors spanning the null space
  if ( n>r ) {
    Z[nopiv,] <- diag(1, n-r, n-r)
    if ( r>0 ) {
      Z[pivcol,] <- -R[1:r,nopiv]
    }
  }
  return (Z) 
  
}



#' Transform matrix A into reduced row echelon form 
#' this function is written along the lines of the rref-matlab function.
#' @param A matrix for which the reduced row echelon form is searched
#' @param tol tolerance to find pivots
#' @param verbose logical, print verbose information
#' @param fractions logical, not used right now. 
#' @return a list of two entries is returned; ret[[1]] is the reduced row echelon form of A, ret[[2]] is the index of columns in which a pivot was found
#' 
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
rref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE){
  ## Written by John Fox
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  m <- nrow(A)
  n <- ncol(A)
  
  i <- 1 # row index
  j <- 1 # column index
  pivcol <- c() # vector of columns in which nozero pivots are found
  while ((i <= m) & (j <= n)){
    # find pivot in column j
    which <- which.max(abs(A[i:m,j])) # column in which pivot is
    k <- i+which-1 #row index, in which pivot is
    pivot <- A[k, j] # pivot of column j
    
    if ( abs(pivot) <= tol ) {
      A[i:m,j] =matrix(0,m-i+1,1) # column is negligible, zero it out
      j <- j+1
    } else {
      # remember column index
      pivcol <- cbind(pivcol,j)
      
      # swap i-th and k-th column
      A[cbind(i,k),j:n] = A[cbind(k, i),j:n];
      
      # divide pivot row by pivot element.
      A[i,j:n] = A[i,j:n]/A[i,j];
      
      # subtract multiples of pivot row from all other rows.
      otherRows <- 1:m
      otherRows <- otherRows[-i]
      for (u in otherRows) {
        A[u,j:n] = A[u,j:n] - A[u,j]*A[i,j:n];
      }
      i = i + 1;
      j = j + 1;
    }
  }
  return (list(A,pivcol))
}

