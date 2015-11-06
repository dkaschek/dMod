## Todo:  Decide on list of attributes that each prdfn should carry



## Methods for class prdlist ------------------------------------------------



#' @export
#' @rdname prdlist
as.prdlist <- function(x, ...) {
  UseMethod("as.prdlist", x)
}

#' @export
#' @rdname prdlist
as.prdlist.list <- function(mylist = NULL, mynames = names(mylist)) {

  if (is.null(mylist)) mylist <- list()

  if (length(mynames) != length(mylist)) stop("names argument has wrong length")

  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- c("prdlist", "list")

  return(mylist)

}


#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  as.prdlist(newlist)
  
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

## Methods for class prdframe ----------------------------
#' @export
#' @rdname plotCombined
plot.prdframe <- function(prediction, data = NULL, ..., scales = "free", facet = "wrap") {
  
  prediction <- list("1" = prediction)
  if(!is.null(data) && is.data.frame(data))
    data <- list("1" = data)
  
  
  plotCombined(prediction = prediction, data = data, ..., scales = scales, facet = facet)
  
}



## Methods for class prdfn ----------------------------------

# print.prdfn <- function(x) {
#   
#   myargs <- args(x)
#   returns <- tail(body(x), 1)
#   
#   no.interest <- c("srcref", "class")
#   attribs <- attributes(x)
#   attribs <- attribs[setdiff(names(attribs), no.interest)]
#   
#   cat("Prediction function:\n")
#   print(myargs)
#   cat("\nFunction returns:\n")
#   print(returns)
#   cat("\nAttributes of the prediction function:\n")
#   for(n in names(attribs)) {
#     cat("\t", n, "\n")
#     cat("\t", attribs[[n]], "\n")
#   }
#   cat("\n")
#   
# }
