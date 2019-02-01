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
    if(any(class(vec1.inner) %in% c("eqnvec", "data.frame"))) {
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
      i <- sapply(l, is.factor)
      l[i] <- lapply(l[i], as.character)
      present.list <- as.list(l)
      missing.names <- setdiff(mynames, names(present.list))
      missing.list <- structure(as.list(rep(NA, length(missing.names))), names = missing.names)
      combined.data <- do.call(function(...) cbind.data.frame(..., stringsAsFactors = FALSE), c(present.list, missing.list))
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
wide2long <- function(out, keep = 1, na.rm = FALSE) {
  
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
  
  outlong <- do.call(rbind, lapply(1:max(c(length(conditions), 1)), function(cond) {
    
    cbind(wide2long.matrix(out[[cond]]), condition = conditions[cond])
    
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
    if (nrow(myout) > 0)
      myout[["condition"]] <- numconditions[cond]
    else
      myout[["condition"]] <- character(0)
    
    return(myout)
    
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



#' Compile one or more prdfn, obsfn or parfn objects
#' 
#' @param ... Objects of class parfn, obsfn or prdfn
#' @param output Optional character of the file to be produced. If several objects were
#' passed, the different C files are all compiled into one shared object file.
#' @param args Additional arguments for the R CMD SHLIB call, e.g. \code{-leinspline}.
#' @param verbose Print compiler output to R command line.
#' @param cores Number of cores used for compilation when several files are compiled.
#' 
#' @export
compile <- function(..., output = NULL, args = NULL, cores = 1, verbose = F) {
  
  objects <- list(...)
  obj.names <- as.character(substitute(list(...)))[-1]
  # Get full list of .c and .cpp files for the obsfn, parfn and prdfn objects in ...
  files <- NULL
  for (i in 1:length(objects)) {
    
    if (inherits(objects[[i]], c("obsfn", "parfn", "prdfn"))) {
      # Get and reset modelname
      filename <- modelname(objects[[i]])
      # Expand modelname by possible endings and check if file exists
      filename <- outer(filename, c("", "_deriv", "_s", "_sdcv", "_dfdx", "_dfdp"), paste0)
      files.obj <- c(paste0(filename, ".c"), paste0(filename, ".cpp"))
      files.obj <- files.obj[file.exists(files.obj)]
      files <- union(files, files.obj)
    }
  }
  
  
  roots <- sapply(files, function(f) {
    l <- strsplit(f, split = ".", fixed = TRUE)[[1]]
    paste(l[1:(length(l)-1)], collapse = ".")
  })
  
  .so <- .Platform$dynlib.ext
  #print(files)
  
  # Sanitize cores on windows
  if (Sys.info()[['sysname']] == "Windows") cores <- 1
  
  #return(files)
  if (is.null(output)) {
    compilation_out <- mclapply(1:length(files), function(i) {
      try(dyn.unload(paste0(roots[i], .so)), silent = TRUE)
      system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", files[i], " ", args), intern = !verbose)
    }, mc.cores = cores, mc.silent = FALSE)
    for (r in roots) dyn.load(paste0(r, .so))
  } else {
    for (i in 1:length(objects)) {
      eval(parse(text = paste0("modelname(", obj.names[i], ") <<- '", output, "'")))
    }
    for (r in roots) try(dyn.unload(paste0(r, .so)), silent = TRUE)
    try(dyn.unload(output), silent = TRUE)
    system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", paste(files, collapse = " "), " -o ", output, .so, " ", args), intern = !verbose)
    dyn.load(paste0(output, .so))
  }
}



#' Determine loaded DLLs available in working directory
#' 
#' @return Character vector with the names of the loaded DLLs available in the working directory
#' @export
getLocalDLLs <- function() {
  
  all.dlls <- getLoadedDLLs()
  is.local <- sapply(all.dlls, function(x) grepl(getwd(), unclass(x)$path, fixed = TRUE))
  names(is.local)[is.local]
  
}


#' Load shared object for a dMod object
#' 
#' Usually when restarting the R session, although all objects are saved in
#' the workspace, the dynamic libraries are not linked any more. \code{loadDLL}
#' is a wrapper for \code{dyn.load} that uses the "modelname" attribute of
#' dMod objects like prediction functions, observation functions, etc. to
#' load the corresponding shared object.
#' 
#' @param ... objects of class prdfn, obsfn, parfn, objfn, ...
#' 
#' @export
loadDLL <- function(...) {
  
  .so <- .Platform$dynlib.ext
  models <- modelname(...)
  files <- paste0(outer(models, c("", "_s", "_sdcv", "_deriv"), paste0), .so)
  files <- files[file.exists(files)]
  
  for (f in files) {
    try(dyn.unload(f), silent = TRUE)
    dyn.load(f)
  }
  
  
  message("The following local files were dynamically loaded: ", paste(files, collapse = ", "))
  
  
  
}


sanitizeCores <- function(cores)  {
  
  max.cores <- parallel::detectCores()
  min(max.cores, cores)
 #  
 # if (Sys.info()[['sysname']] == "Windows") cores <- 1
 # return(cores)
  
}

sanitizeConditions <- function(conditions) {
  
  new <- str_replace_all(conditions, "[[:punct:]]", "_")
  new <- str_replace_all(new, "\\s+", "_")
  return(new)
  
}

sanitizePars <- function(pars = NULL, fixed = NULL) {
  
  # Convert fixed to named numeric
  if (!is.null(fixed)) fixed <- structure(as.numeric(fixed), names = names(fixed))
  
  # Convert pars to named numeric
  if (!is.null(pars)) {
    pars <- structure(as.numeric(pars), names = names(pars))
    # remove fixed from pars
    pars <- pars[setdiff(names(pars), names(fixed))]
  }
    
  
  return(list(pars = pars, fixed = fixed))
  
}


sanitizeData <- function(x, required = c("name", "time", "value"), imputed = c(sigma = NA, lloq = -Inf)) {
  
  all.names <- names(x)
  
  missing.required <- setdiff(required, all.names)
  missing.imputed <- setdiff(names(imputed), all.names)
  
  if (length(missing.required) > 0)
      stop("These mandatory columns are missing: ", paste(missing.required, collapse = ", "))
  
  if (length(missing.imputed) > 0) {
    
    for (n in missing.imputed) x[[n]] <- imputed[n]
    
  }
  
  list(data = x, columns = c(required, names(imputed)))
  
}


#' Print list of dMod objects in .GlobalEnv
#' 
#' @description Lists the objects for a set of classes.
#'   
#' @param classlist List of object classes to print.
#' @param envir Alternative environment to search for objects.
#' @examples 
#' \dontrun{
#' lsdMod()
#' lsdMod(classlist = "prdfn", envir = environment(obj)) 
#' }
#' 
#' @export
lsdMod <- function(classlist = c("odemodel", "parfn", "prdfn", "obsfn", "objfn", "datalist"), envir = .GlobalEnv){
  glist <- as.list(envir)
  out <- list()
  for (a in classlist) {
    flist <- which(sapply(glist, function(f) any(class(f) == a)))
    out[[a]] <- names(glist[flist])
    #cat(a,": ")
    #cat(paste(out[[a]], collapse = ", "),"\n")
  }
  
  unlist(out)
  
  
  
}





#' Select attributes.
#' 
#' @description Select or discard attributes from an object.
#'   
#' @param x The object to work on
#' @param atr An optional list of attributes which are either kept or removed. 
#'   This parameter defaults to dim, dimnames, names,  col.names, and row.names.
#' @param keep For keep = TRUE, atr is a positive list on attributes which are 
#'   kept, for keep = FALSE, \option{atr} are removed.
#'   
#' @return x with selected attributes.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
attrs <- function(x, atr = NULL, keep = TRUE) {

  if (is.null(atr)) {
    atr <- c("class", "dim", "dimnames", "names", "col.names", "row.names")
  }
  
  xattr <- names(attributes(x))
  if (keep == TRUE) {
    attributes(x)[!xattr %in% atr] <- NULL
  } else {
    attributes(x)[xattr %in% atr] <- NULL
  }
  
  return(x)
}




#' Print object and its "default" attributes only.
#' 
#' @param x Object to be printed
#' @param list_attributes Prints the names of all attribute of x, defaults to 
#'   TRUE
#'   
#' @details Before the \option{x} is printed by print.default, all its arguments
#'   not in the default list of \code{\link{attrs}} are removed.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, 
#'   \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
print0 <- function(x, list_attributes = TRUE ) {
  if (list_attributes == TRUE) {
    cat("List of all attributes: ", names(attributes(x)), "\n")
  }
  
  print.default(attrs(x))
}


