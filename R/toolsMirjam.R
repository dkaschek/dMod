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


#' Named repititions
#' 
#' @description Wrapper on rep() to input names instead of length.
#'   
#' @param x Value to be repeated.
#' @param names List of names.
#'   
#' @export
repWithNames <- function(x, names){
  repnum <- rep(x,length(names))
  names(repnum) <- names
  return(repnum)
}
