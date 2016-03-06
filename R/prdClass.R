
## Methods for class prdlist ------------------------------------------------



#' @export
#' @rdname prdlist
as.prdlist <- function(x, ...) {
  UseMethod("as.prdlist", x)
}

#' @export
#' @param x list of prediction frames
#' @param names character vector, the list names, e.g. the names of the experimental
#' @rdname prdlist
as.prdlist.list <- function(x = NULL, names = NULL) {

  if (is.null(x)) x <- list()
  if (is.null(names)) mynames <- names(x) else mynames <- names 

  # if (length(mynames) != length(x)) stop("names argument has wrong length")

  ## Prepare output
  names(x) <- mynames
  class(x) <- c("prdlist", "list")

  return(x)

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
plot.prdlist <- function(prediction, data = NULL, ..., scales = "free", facet = "wrap", transform = NULL) {
  
  if (is.null(names(prediction))) names(prediction) <- 1:length(prediction)
  if (!is.null(data) && is.null(names(data))) names(data) <- 1:length(data)
  
  plotCombined(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}

## Methods for class prdframe ----------------------------
#' @export
#' @rdname plotCombined
plot.prdframe <- function(prediction, data = NULL, ..., scales = "free", facet = "wrap", transform = NULL) {
  
  prediction <- list("1" = prediction)
  if(!is.null(data) && is.data.frame(data))
    data <- list("1" = data)
  
  
  plotCombined(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}



## Methods for class prdfn ----------------------------------

#' @export
print.prdfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Prediction function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.prdfn <- function(x, ...) {
  
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        events = attr(mappings[[C]], "events"),
        forcings = attr(mappings[[C]], "forcings"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    print(output, ...)
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}

#' @export
print.obsfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Observation function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.obsn <- function(x, ...) {
  
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        states = attr(mappings[[C]], "states"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    print(output, ...)
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}


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
