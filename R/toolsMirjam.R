#' Print without attributes
#' 
#' @param x The object to be printed out
#' @param list_attributes prints a list of attribute names if TRUE (default=TRUE)
#' @details To suppress the printout of attributes like "deriv". 
#' @export
print0 <- function(x, list_attributes=TRUE ) {
  attributes_all <- names(attributes(x))
  attributes_rm <- attributes_all[!(attributes_all %in% c("dim","names","dimnames","row.names","col.names"))]
  attributes(x)[attributes_rm] <- NULL
  print.default(x)
  if(list_attributes)
    cat("Attributes:",attributes_all)
}

#' Find common elements in lists
#' @description generalisation of \link{intersect} to a list of vectors.
#' @param list contains a set of vectors.
#' @param byNames if set TRUE common names are checked instead of common values (default: TRUE) 
#' @examples testList <-list(c(a=1,b=5,c=3,e=5), c(d=3,b=1,q=1,c=5,i=2)) 
#' intersectList(testList,FALSE) 
#' intersectList(testList,TRUE)
intersectList <- function(list, byNames=TRUE){
  inter <- list[[1]]
  if(byNames)
    inter <- names(inter)
  lapply(list[-1], function(l){
    if(byNames)
      l <- names(l)
    inter <<- intersect(inter,l)
  })
  return(inter)
}

