## Todo:  Decide on list of attributes that each prdfn should carry



## Methods for class prdlist ------------------------------------------------

#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  prdlist(newlist)
  
}

#' @export
"[.prdlist" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("prdlist", "list")
  return(out)
}

#' @export
#' @rdname plotCombined
plot.prdlist <- function(prediction, data = NULL, ..., scales = "free", facet = "wrap") {
  
  plotCombined(prediction = prediction, data = data, ..., scales = scales, facet = facet)
  
}