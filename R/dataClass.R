
## Class "datalist" and its constructor ------------------------------------------


#' @param x object of class \code{data.frame} or \code{list} 
#' @return Object of class \link{datalist}
#' @export
#' @example inst/examples/datalist.R 
#' @rdname datalist
as.datalist <- function(x, ...) {
  UseMethod("as.datalist", x)
}

#' @export
#' @param split.by vector of columns names which yield a unique identifier (conditions). If NULL, all
#' columns except for the expected standard columns "name", "time", "value" and "sigma" will be
#' selected.
#' @rdname datalist
as.datalist.data.frame <- function(x, split.by = NULL, ...) {
  
  dataframe <- x
  
  #remaining.names <- setdiff(names(dataframe), split.by)
  all.names <- colnames(dataframe)
  standard.names <- c("name", "time", "value", "sigma")
  if (is.null(split.by)) split.by <- setdiff(all.names, standard.names)
  
  
 
  conditions <- lapply(split.by, function(n) dataframe[, n])
  splits <- do.call(paste, c(conditions, list(sep = "_")))
  
  # condition grid
  conditionframe <- dataframe[!duplicated(splits), split.by, drop = FALSE]
  rownames(conditionframe) <- splits[!duplicated(splits)]
  
  
  # data list output
  dataframe <- cbind(data.frame(condition = splits), dataframe[, standard.names])
  out <- lapply(unique(splits), function(s) dataframe[dataframe[, 1] == s, -1])
  
  names(out) <- as.character(unique(splits))
  
  out <- as.datalist(out)
  attr(out, "condition.grid") <- conditionframe
  return(out)
  
}

#' @export
#' @param names optional names vector, otherwise names are taken from \code{mylist}
#' @rdname datalist
as.datalist.list <- function(x, names = NULL, ...) {
  
  mylist <- x
  
  ## Check properties
  if (is.null(names)) mynames <- names(mylist) else mynames <- names
  is.data.frame <- sapply(mylist, class) == "data.frame"
  if (!all(is.data.frame)) stop("list of data.frame expected")

  correct.names <- c("name", "time", "value", "sigma")
  have.correct.names <- sapply(mylist, function(d) all(correct.names %in% colnames(d)))
  if (all(have.correct.names)) {
    mylist <- lapply(mylist, function(d) d[, correct.names])
  } else {
    stop(paste("data.frames should have names:", correct.names, collapse = " "))    
  }


  if (length(mynames) != length(mylist)) stop("names argument has wrong length")

  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- c("datalist", "list")
  attr(mylist, "condition.grid") <- data.frame(condition = mynames, row.names = mynames)

  return(mylist)

}


## Methods for class datalist ---------------------------------------

#' @export
print.datalist <- function(x, ...) {
  datalist <- x
  for(n in names(datalist)) {
    cat(n, ":\n", sep = "")
    print(datalist[[n]])
  }
}

# Subset of all datalist entries
#' @export
subset.datalist <- function(x, ...){
  datalist <- lapply(x, function(i) subset(i, ...)) 
  return(as.datalist(datalist))
}

#' @export
"[.datalist" <- function(x, ...) {
  condition.grid <- attr(x, "condition.grid")
  out <- unclass(x)[...]
  attr(out, "condition.grid") <- condition.grid
  n <- names(out)
  if (!is.null(n)) {
    attr(out, "condition.grid") <- condition.grid[n, , drop = FALSE]
  }
  class(out) <- c("datalist", "list")
  return(out)
}

#' Plot a list data points
#' 
#' @param x Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{subset}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @export
plot.datalist <- function(x, ..., scales = "free", facet = "wrap") {
  
  data <- x
  if (is.null(names(data))) names(data) <- paste0("C", 1:length(data))
  plotCombined(prediction = NULL, data = data, ..., scales = scales, facet = facet)
  
}

#' Coerce to a Data Frame
#' 
#' @param x any R object
#' @return a data frame
#' @rdname as.data.frame.dMod
#' @export
as.data.frame.datalist <- function(x, ...) {
  
  data <- x
  condition.grid <- attr(x, "condition.grid")
  
  data <- lbind(data)
  if (!is.null(condition.grid)) {
    for (C in colnames(condition.grid)) {
      data[, C] <- condition.grid[as.character(data$condition), C]
    }
  }
  
  return(data)
  
}