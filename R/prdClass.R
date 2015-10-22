## Todo:  Decide on list of attributes that each prdfn should carry



## Methods for class prdlist ------------------------------------------------

#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  prdlist(newlist)
  
}
