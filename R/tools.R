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
    
    #myout <- out[[cond]]
    #timename <- colnames(myout)[1]
    #allnames <- colnames(myout)[-1]
    #times <- myout[,1]
    #values <- unlist(myout[,allnames])
    #myoutlong <- data.frame(time = times, 
    #                        name = rep(allnames, each=length(times)), 
    #                        value = as.numeric(values), 
    #                        condition = cond)
    #colnames(myoutlong)[1] <- timename
    #return(myoutlong)
    
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





#' Evaluation of algebraic expressions defined by characters
#' 
#' @param x Name character vector, the algebraic expressions
#' @param compile Logical. The function is either translated into a C file to be compiled or is
#' evaluated in raw R.
#' @return A prediction function \code{f(mylist, attach.input = FALSE)} where \code{mylist} is a list of numeric 
#' vectors that can
#' be coerced into a matrix. The names correspond to the symbols used in the algebraic expressions. 
#' The argument \code{attach.input} determines whether \code{mylist} is attached to the output.
#' The function \code{f} returns a matrix.
#' @examples 
#' \dontrun{
#' myfun <- funC0(c(y = "a*x^4 + b*x^2 + c"))
#' out <- myfun(list(a = -1, b = 2, c = 3, x = seq(-2, 2, .1)), attach.input = TRUE)
#' plot(out[, "x"], out[, "y"])
#' }
#' 
#' @export
funC0 <- function(x, compile = FALSE, modelname = NULL) {
    
  # Get symbols to be substituted by x[] and y[]
  outnames <- names(x)
  innames <- getSymbols(x)
  
  x.new <- paste0(x, collapse = ", ")
  x.new <- paste0("list(", x.new, ")")
  x.expr <- parse(text = x.new)
  
  ## Compiled version based on inline package
  ## Non-compiled version based on with() and eval()
  if(compile) {
    
    # Do the replacement to obtain C syntax
    x <- replaceOperation("^", "pow", x)
    x <- replaceSymbols(innames, paste0("x[", (1:length(innames))-1, "+i* *k]"), x)
    names(x) <- paste0("y[", (1:length(outnames)) - 1, "+i* *l]")
    
    # Paste into equation
    x <- x[x != "0"]
    expr <- paste(names(x), "=", x, ";")
    
    # Put equation into C function
    if(is.null(modelname)) {
      funcname <- paste0("funC0_", paste(sample(c(0:9, letters), 8, replace = TRUE), collapse = ""))
    } else {
      funcname <- modelname
    }
    body <- paste(
      "#include <R.h>\n", 
      "#include <math.h>\n", 
      "void", funcname, "( double * x, double * y, int * n, int * k, int * l ) {\n",
      "for(int i = 0; i< *n; i++) {\n",
      paste(expr, collapse="\n"),
      "\n}\n}"
    )
    
    filename <- paste(funcname, "c", sep = ".")
    sink(file = filename)
    cat(body)
    sink()
    system(paste0(R.home(component="bin"), "/R CMD SHLIB ", filename))
    .so <- .Platform$dynlib.ext
    dyn.load(paste0(funcname, .so))
    
    # Generate the C function by the inline package
    #myCfun <- inline::cfunction(sig=c(x = "double", y = "double", n = "integer", k = "integer", l = "integer"),
    #                            body=body,
    #                            language="C",
    #                            convention=".C"
    #)
    
    
    # Generate output function
    myRfun <- function(x, attach.input = FALSE) {
      
      # Translate the list into matrix and then into vector
      M <- do.call(rbind, x[innames])
      if(length(M) == 0) M <- matrix(0)
      x <- as.double(as.vector(M))
      
      # Get integers for the array sizes
      n <- as.integer(dim(M)[2])
      k <- as.integer(length(innames))
      if(length(k) == 0) k <- as.integer(0)
      l <- as.integer(length(outnames))
      
      
      # Initialize output vector
      y <- double(l*n)
      
      # Evaluate C function and write into matrix
      loadDLL(func = funcname, cfunction = funcname)
      out <- matrix(.C(funcname, x = x, y = y, n = n, k = k, l = l)$y, nrow=length(outnames), ncol=n)
      rownames(out) <- outnames
      
      rownames(M) <- innames
      if(attach.input)
        out <- rbind(M, out)
        
      
      return(t(out))    
      
    }
    
    
    
  } else {
    
    # Generate output function
    myRfun <- function(x, attach.input = FALSE) {
      
      # Translate the list into matrix and then into vector
      M <- do.call(rbind, x[innames])
      if(length(M) == 0) M <- matrix(0)
      
      out.list <- with(x, eval(x.expr))
      out.matrix <- do.call(cbind, out.list)
      colnames(out.matrix) <- outnames
      rownames(out.matrix) <- NULL
      
      rownames(M) <- innames
      if(attach.input)
        out.matrix <- cbind(t(M), out.matrix)
      
      return(out.matrix)
      
    }
    
  }
  
  
  attr(myRfun, "equations") <- x
  
  return(myRfun)
  
}

