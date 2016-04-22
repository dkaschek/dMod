#' Compare two objects and return differences
#' 
#' Works eigher on a list or on two arguments. In case of a list,
#' comparison is done with respect to a reference entry. Besides the
#' objects themselves also some of their attributes are compared,
#' i.e. "equations", "parameters" and "events" and "forcings".
#' 
#' @param vec1 object of class \link{eqnvec}, \code{character} or
#' \code{data.frame}. Alternatively, a list of such objects.
#' @param vec2 same as vec1. Not used if vec1 is a list.
#' @param reference numeric of length one, the reference entry.
#' @param ... arguments going to the corresponding methods
#' @return \code{data.frame} or list of data.frames with the differences. 
#' 
#' @export
#' @examples
#' ## Compare equation vectors
#' eq1 <- eqnvec(a = "-k1*a + k2*b", b = "k2*a - k2*b")
#' eq2 <- eqnvec(a = "-k1*a", b = "k2*a - k2*b", c = "k2*b")
#' compare(eq1, eq2)
#' 
#' ## Compare character vectors
#' c1 <- c("a", "b")
#' c2 <- c("b", "c")
#' compare(c1, c2)
#' 
#' ## Compare data.frames
#' d1 <- data.frame(var = "a", time = 1, value = 1:3, method = "replace")
#' d2 <- data.frame(var = "a", time = 1, value = 2:4, method = "replace")
#' compare(d1, d2)
#' 
#' ## Compare structures like prediction functions
#' fn1 <- function(x) x^2
#' attr(fn1, "equations") <- eq1
#' attr(fn1, "parameters") <- c1
#' attr(fn1, "events") <- d1
#' 
#' fn2 <- function(x) x^3
#' attr(fn2, "equations") <- eq2
#' attr(fn2, "parameters") <- c2
#' attr(fn2, "events") <- d2
#' 
#' mylist <- list(f1 = fn1, f2 = fn2)
#' compare(mylist)

compare <- function(vec1, ...) {
  UseMethod("compare", vec1)
}

#' @export
#' @rdname compare
compare.list <- function(vec1, vec2 = NULL, reference = 1, ...) {
  
  index <- (1:length(vec1))[-reference]
  diffable.attributes <- c("equations", "parameters", "forcings", "events")
  
  
  out.total <- lapply(index, function(i) {
    
    # Compare objects if possible
    vec1.inner <- vec1[[reference]]
    vec2.inner <- vec1[[i]]
    out1 <- NULL
    if(class(vec1.inner) %in% c("eqnvec", "data.frame")) {
      out1 <- list(compare(vec1.inner, vec2.inner))
      names(out1) <- "object"
    }
      
    # Compare comparable attributes of the object if available
    out2 <- NULL
    attributes1 <- attributes(vec1.inner)[diffable.attributes]
    attributes2 <- attributes(vec2.inner)[diffable.attributes]
    slots <- names(attributes1)[!is.na(names(attributes1))]
    out2 <- lapply(slots, function(n) {
      compare(attributes1[[n]], attributes2[[n]])
    })
    names(out2) <- slots
    
    c(out1, out2)
    
  })
  names(out.total) <- names(vec1)[index]
  
  ## Do resorting of the list
  innernames <- names(out.total[[1]])
  out.total <- lapply(innernames, function(n) {
    out <- lapply(out.total, function(out) out[[n]])
    out[!sapply(out, is.null)]
  })
  names(out.total) <- innernames
  
  
  return(out.total)
  
  
  
}

#' @export
#' @rdname compare
compare.character <- function(vec1, vec2 = NULL, ...) {
  missing <- setdiff(vec1, vec2)
  additional <- setdiff(vec2, vec1)
  
  out <- do.call(rbind, 
          list(different = NULL, 
               missing = data.frame(name = missing), 
               additional = data.frame(name = additional)
          )
  )
  
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
  
}

#' @export
#' @rdname compare
compare.eqnvec <- function(vec1, vec2 = NULL, ...) {

  names1 <- names(vec1)
  names2 <- names(vec2)
  
  missing <- setdiff(names1, names2)
  additional <- setdiff(names2, names1)
  joint <- intersect(names1, names2)
  
  # Compare joint equations
  v1 <- format(vec1)
  v2 <- format(vec2)
  not.coincide <- which(as.character(v1[joint]) != as.character(v2[joint]))
  
  different <- data.frame(name = names(v2[not.coincide]), equation = as.character(v2[not.coincide]))
  missing <- data.frame(name = names(v2[missing]), equation = as.character(v2[missing]))
  additional <- data.frame(name = names(v2[additional]), equation = as.character(v2[additional]))
  
  out <- do.call(rbind, list(different = different, missing = missing, additional = additional))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}

#' @export
#' @rdname compare
compare.data.frame <- function(vec1, vec2 = NULL, ...) {
  
  additional <- !duplicated(rbind(vec1, vec2))[-(1:nrow(vec1))]
  missing <- !duplicated(rbind(vec2, vec1))[-(1:nrow(vec2))]
  
  out <- do.call(rbind, list(different = character(0), missing = vec1[missing, ], additional = vec2[additional, ]))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}


#' Combine several data.frames by rowbind
#' 
#' @param ... data.frames or matrices with not necessarily overlapping colnames
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
#' @export
combine <- function(...) {
  
  # List of input data.frames
  mylist <- list(...)
  # Remove empty slots
  is.empty <- sapply(mylist, is.null)
  mylist <- mylist[!is.empty]
  
  mynames <- unique(unlist(lapply(mylist, function(S) colnames(S))))
  
  mylist <- lapply(mylist, function(l) {
    
    if(is.data.frame(l)) {
      present.list <- as.list(l)
      missing.names <- setdiff(mynames, names(present.list))
      missing.list <- structure(as.list(rep(NA, length(missing.names))), names = missing.names)
      combined.data <- do.call(cbind.data.frame, c(present.list, missing.list))
      rownames(combined.data) <- rownames(l)
    }
    if(is.matrix(l)) {
      present.matrix <- as.matrix(l)
      missing.names <- setdiff(mynames, colnames(present.matrix))
      missing.matrix <- matrix(0, nrow = nrow(present.matrix), ncol = length(missing.names), 
                             dimnames = list(NULL, missing.names))
      combined.data <- submatrix(cbind(present.matrix, missing.matrix), cols = mynames)
      rownames(combined.data) <- rownames(l)
    }
    
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
  
 M[rows, cols, drop = FALSE] 
  
  # myrows <- (structure(1:nrow(M), names = rownames(M)))[rows]
  # mycols <- (structure(1:ncol(M), names = colnames(M)))[cols]
  # 
  # if(any(is.na(myrows)) | any(is.na(mycols))) stop("subscript out of bounds")
  # 
  # matrix(M[myrows, mycols], 
  #        nrow = length(myrows), ncol = length(mycols), 
  #        dimnames = list(rownames(M)[myrows], colnames(M)[mycols]))

}


#' Embed two matrices into one blockdiagonal matrix
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character containing M and N as upper left and lower right block
#' @examples
#' M <- matrix(1:9, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
#' N <- matrix(1:4, 2, 2, dimnames = list(LETTERS[1:2], LETTERS[1:2]))
#' blockdiagSymb(M, N)
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
  colnames(result) <- c(colnames(M), colnames(N))
  rownames(result) <- c(rownames(M), rownames(N))
  
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
  if (any(duplicated(allnames))) warning("Found duplicated colnames in out. Duplicates were removed.")
  times <- out[,keep]
  ntimes <- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times, 
                        name = factor(rep(allnames, each = ntimes), levels = allnames), 
                        value = as.numeric(values))
  colnames(outlong)[1:length(keep)] <- timenames
  
  if (na.rm) outlong <- outlong[!is.na(outlong$value),]
  
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
  #numconditions <- suppressWarnings(as.numeric(conditions))
  #
  # if(!any(is.na(numconditions))) 
  #   numconditions <- as.numeric(numconditions) 
  # else 
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





