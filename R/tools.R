#' Combine several data.frames by rowbind
#' 
#' @param ... data.frames with not necessarily overlapping colnames
#' @details This function is useful when separating models into independent csv model files,
#' e.g.~a receptor model and several downstream pathways. Then, the models can be recombined 
#' into one model by \code{combine()}.
#' 
#' @return A \code{data.frame}
#' @export
#' @examples
#' data1 <- data.frame(Description = "reaction 1", Rate = "k1*A", A = -1, B = 1)
#' data2 <- data.frame(Description = "reaction 2", Rate = "k2*B", B = -1, C = 1)
#' combine(data1, data2)
combine <- function(...) {
  
  # List of input data.frames
  mylist <- list(...)
  # Remove empty slots
  is.empty <- sapply(mylist, is.null)
  mylist <- mylist[!is.empty]
  
  mynames <- unique(unlist(lapply(mylist, function(S) colnames(S))))
  
  mylist <- lapply(mylist, function(l) {
    present.list <- as.list(l)
    missing.names <- setdiff(mynames, names(present.list))
    missing.list <- structure(as.list(rep(NA, length(missing.names))), names = missing.names)
    combined.data <- do.call(cbind.data.frame, c(present.list, missing.list))
    return(combined.data)
  })
  
  out <- do.call(rbind, mylist)
  
  return(out)
  
  
}

#' Submatrix of a matrix returning ALWAYS a matrix
#' 
#' @param M matrix
#' @param rows Index vector
#' @param cols Index vector
#' @return The matrix \code{M[rows, cols]}, keeping/adjusting attributes like ncol nrow and dimnames.
#' @export
submatrix <- function(M, rows = 1:nrow(M), cols = 1:ncol(M)) {
  
  
  
  myrows <- (structure(1:nrow(M), names = rownames(M)))[rows]
  mycols <- (structure(1:ncol(M), names = colnames(M)))[cols]
  
  if(any(is.na(myrows)) | any(is.na(mycols))) stop("subscript out of bounds")
  
  matrix(M[myrows, mycols], 
         nrow = length(myrows), ncol = length(mycols), 
         dimnames = list(rownames(M)[myrows], colnames(M)[mycols]))
}


#' Embed two matrices into one blockdiagonal matrix
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character containing M and N as upper left and lower right block
#' @export
blockdiagSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  A <- matrix(0, ncol=dim(N)[2], nrow=dim(M)[1])
  B <- matrix(0, ncol=dim(M)[2], nrow=dim(N)[1])
  result <- rbind(cbind(M, A), cbind(B, N))
  return(result)
  
}



#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export
wide2long <- function(out, keep, na.rm) {
  
  UseMethod("wide2long", out)
  
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.data.frame
#' @export
wide2long.data.frame <- function(out, keep = 1, na.rm = FALSE) {
  
  wide2long.matrix(out, keep = keep, na.rm = na.rm)
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.matrix
#' @export
wide2long.matrix <- function(out, keep = 1, na.rm = FALSE) {
  
  timenames <- colnames(out)[keep]
  allnames <- colnames(out)[-keep]
  times <- out[,keep]
  ntimes<- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times, name = rep(allnames, each=ntimes), value = as.numeric(values))
  colnames(outlong)[1:length(keep)] <- timenames
  
  if(na.rm) outlong <- outlong[!is.na(outlong$value),]
  
  return(outlong)
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.list
#' @export
wide2long.list <- function(out, keep = 1, na.rm = FALSE) {
  
  conditions <- names(out)
  numconditions <- suppressWarnings(as.numeric(conditions))
  
  if(!any(is.na(numconditions))) 
    numconditions <- as.numeric(numconditions) 
  else 
    numconditions <- conditions
  
  
  outlong <- do.call(rbind, lapply(1:length(conditions), function(cond) {
    
    cbind(wide2long.matrix(out[[cond]]), condition = numconditions[cond])
    
  }))
  
  
  
  return(outlong)
  
}


#' Translate long to wide format (inverse of wide2long.matrix) 
#' 
#' @param out data.frame in long format 
#' @return data.frame in wide format 
#' @export
long2wide <- function(out) {
  
  timename <- colnames(out)[1]
  times <- unique(out[,1])
  allnames <- unique(as.character(out[,2]))
  M <- matrix(out[,3], nrow=length(times), ncol=length(allnames))
  M <- cbind(times, M)
  colnames(M) <- c(timename, allnames)
  
  return(M)
  
}


#' Bind named list of data.frames into one data.frame
#' 
#' @param mylist A named list of data.frame. The data.frames are expected to have the same structure.
#' @details Each data.frame ist augented by a "condition" column containing the name attributed of
#' the list entry. Subsequently, the augmented data.frames are bound together by \code{rbind}.
#' @return data.frame with the originial columns augmented by a "condition" column.
#' @export
lbind <- function(mylist) {
  
  conditions <- names(mylist)
  numconditions <- suppressWarnings(as.numeric(conditions))
  
  if(!any(is.na(numconditions))) 
    numconditions <- as.numeric(numconditions) 
  else 
    numconditions <- conditions

  
  outlong <- do.call(rbind, lapply(1:length(conditions), function(cond) {
    
    myout <- mylist[[cond]]
    myoutlong <- cbind(myout, condition = numconditions[cond])
    
    return(myoutlong)
    
  }))
  
  return(outlong)
  
}

#' Alternative version of expand.grid
#' @param seq1 Vector, numeric or character
#' @param seq2 Vector, numeric or character
#' @return Matrix ob combinations of elemens of \code{seq1} and \code{seq2}
expand.grid.alt <- function(seq1, seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}

#' Load a template file in the editor
#' 
#' @param i Integer, choose a template to be loaded
#' @details Possible templates are:
#' i = 1: Do parameter estimation in a dynamic model with fixed forcings
#' @export
loadTemplate <- function(i = 1) {
  
  path <- path.package("dMod")
  if(i == 1) {
    system(paste0("cp ", path, "/templates/R2CTemplate.R mymodel.R"))
    file.edit("mymodel.R")
  }
  
}





