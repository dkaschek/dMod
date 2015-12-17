#' Plot a parameter list.
#' 
#' @param pl fitlist obtained from mstrust
#' @param path print path of parameters from initials to convergence. For this
#'   option to be TRUE \code{\link{mstrust}} must have had the option
#'   \option{blather}.
#' 
#' @details If path=TRUE:        
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
plot.parlist <- function(pl, path = FALSE) {
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
#' @param verbose option for rref-function below
#' @param fractions option for rref-function below
#' 
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
nullZ <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE) {
  R <- rref(A) # compute reduced row echelon form of A
  n <- ncol(A) # number of columns of A
  r <- Matrix::rankMatrix(R) # rank of reduced row echelon form
  if ( r>0 ) {
    pivcol <- 1:r # columns in which a pivot was found (this is by definition equivalent to the first r columns)
  } else {
    pivcol <- c() # for r=0 there are no pivots
  }
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



#' Transform matrix A into reduced row echelon form (by Gauss-Jordan elimination)  
#' this function is taken from stackoverflow.com/questions/3126759/reduced-row-echelon-form and tested in parallel to matlab results.
#' @param A matrix for which the reduced row echelon form is searched
#' @param tol tolerance to find pivots
#' @param verbose option for printing intermediate steps
#' @param fractions try to express nonintegers as rational numbers
#' 
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
rref <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, fractions=FALSE){
  ## Written by John Fox
  if (fractions) {
    mass <- require(MASS)
    if (!mass) stop("fractions=TRUE needs MASS package")
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  for (i in 1:min(c(m, n))){
    col <- A[,i]
    col[1:n < i] <- 0
    # find maximum pivot in current column at or below current row
    which <- which.max(abs(col))
    pivot <- A[which, i]
    if (abs(pivot) <= tol) next     # check for 0 pivot
    if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows
    A[i,] <- A[i,]/pivot            # pivot
    row <- A[i,]
    A <- A - outer(A[,i], row)      # sweep
    A[i,] <- row                    # restore current row
    if (verbose)
      if (fractions) print(fractions(A))
    else print(round(A,round(abs(log(tol,10)))))
  }
  for (i in 1:n)
    if (max(abs(A[i,1:m])) <= tol)
      A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
  if (fractions) {
    return(fractions (A))
  } else {
    return(round(A, round(abs(log(tol,10)))))
  }
}
