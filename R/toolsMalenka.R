#' Plot a parameter list.
#' 
#' @param x fitlist obtained from mstrust
#' @param ... additional arguments
#' @param path print path of parameters from initials to convergence. For this
#'   option to be TRUE \code{\link{mstrust}} must have had the option
#'   \option{blather}.
#' 
#' @details If path=TRUE:        
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
plot.parlist <- function(x, path = FALSE, ...) {
  
  pl <- x
  
  index <- do.call(rbind, lapply(pl, function(l) l$converged))
  fl <- pl[index]
  if (!path) {
    initPar <- do.call(rbind, lapply(fl, function(l) l$parinit))
    convPar <- do.call(rbind, lapply(fl, function(l) l$argument))
    
    ddata <- data.frame(cbind(matrix(initPar, ncol = 1), matrix(convPar, ncol = 1) ))
    ddata <- cbind(rep(colnames(initPar), each = nrow(initPar)), ddata, 1)
    names(ddata) <- c("parameter","x","y","run")
    
    #plot initial vs converged parameter values
    ggplot(data=ddata)+facet_wrap(~ parameter)+geom_point(aes(x=x,y=y))
  } else {
    if (!any (names(fl[[1]]) == "argpath")){
      stop("No path information in the output of mstrust. Restart mstrust with option blather.")
    }
    parNames <- names(fl[[1]]$parinit)
    
    pathPar <- do.call(rbind, mapply(function(l, idx) {
      mParPath <- as.data.frame(matrix(l$argpath, ncol = 1))
      mParPath <- cbind(rep(parNames,each = nrow(l$argpath), times = 1), rep(1:nrow(l$argpath), length(parNames)), mParPath, as.character(idx))
    }, l = fl, idx = 1:length(fl), SIMPLIFY = FALSE))
    names(pathPar) <- c("parameter", "iteration", "path", "idx")
    ggplot(data=pathPar)+geom_line(aes(x=iteration,y=path,colour=idx))+facet_wrap(~ parameter)
  }
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
