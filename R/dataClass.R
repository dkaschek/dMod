
## Class "datalist" and its constructor ------------------------------------------


#' @param x object of class \code{data.frame} or \code{list}. Data frames are required to
#' provide "name", "time" and "value" as columns. Columns "sigma" and "lloq" can be provided.
#' If "sigma" and "lloq" are missing, they
#' are imputed with \code{NA} and \code{-Inf}, respectively. 
#' @return Object of class \link{datalist}
#' @export
#' @example inst/examples/datalist.R
#' @rdname datalist
as.datalist <- function(x, ...) {
  UseMethod("as.datalist", x)
}

#' @export
#' @param split.by vector of columns names which yield a unique identifier (conditions). If NULL, all
#' columns except for the expected standard columns "name", "time", "value", "sigma" and "lloq" will be
#' selected.
#' @param keep.covariates vector of additional column names which should be kept in the condition.grid.
#' @rdname datalist
as.datalist.data.frame <- function(x, split.by = NULL, keep.covariates = NULL, ...) {

  # Sanitize data and get names
  x <- sanitizeData(x)
  dataframe <- x[["data"]]
  standard.names <- x[["columns"]]
  all.names <- colnames(dataframe)
  
  # Get splitting information
  if (is.null(split.by)) split.by <- setdiff(all.names, standard.names)
  conditions <- lapply(split.by, function(n) dataframe[, n])
  splits <- do.call(paste, c(conditions, list(sep = "_")))


  # condition grid
  conditionframe <- dataframe[!duplicated(splits), union(split.by, keep.covariates), drop = FALSE]
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
#' @param condition.grid Optionally, to manually specify a condition.grid
#' @rdname datalist
as.datalist.list <- function(x, names = NULL, ..., condition.grid = attr(x, "condition.grid")) {

  mylist <- x


  ## Check properties
  if (is.null(names)) mynames <- names(mylist) else mynames <- names
  is.data.frame <- sapply(mylist, class) == "data.frame"
  if (!all(is.data.frame)) stop("list of data.frame expected")

  # Sanitize data in list
  mylist <- lapply(mylist, function(x) sanitizeData(x)[["data"]])

  if (length(mynames) != length(mylist)) stop("names argument has wrong length")

  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- c("datalist", "list")


  if (is.null(condition.grid)) {
    condition.grid <- data.frame(condition = mynames, row.names = mynames)
  }
  attr(mylist, "condition.grid") <- condition.grid

  return(mylist)

}


## Methods for class datalist ---------------------------------------

#' @param value The new condition names of the datalist and its condition.grid
#' @export
#' @rdname datalist
"names<-.datalist" <- function(x, value) {
  x <- unclass(x)
  x <- base::`names<-`(x, value)
  attr(x, "condition.grid") <- base::`rownames<-`(attr(x, "condition.grid"), value)
  return(as.datalist(x))
}

#' @export
#' @rdname datalist
is.datalist <- function(x) {
  inherits(x, "datalist")
}

#' @export
#' @rdname datalist
c.datalist <- function(...) {
  dlist <- lapply(list(...), unclass)
  
  condition.grids <- lapply(dlist, function(i) attr(i, "condition.grid"))
  mycg <- Reduce(dMod::combine, condition.grids)
  
  dlist <- Reduce(c, lapply(dlist, function(i) {`attr<-`(i, "condition.grid", NULL)}))
  attr(dlist, "condition.grid") <-  mycg
  class(dlist) <- "datalist"
  return(dlist)
}


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
#' @param ... Further arguments going to \code{dplyr::filter}.
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
  plotCombined.prdlist(prediction = NULL, data = data, ..., scales = scales, facet = facet)

}

#' @export
#' @rdname plotData
#' @importFrom dplyr filter
plotData.datalist <- function(data, ..., scales = "free", facet = "wrap", transform = NULL) {

  rownames_to_condition <- function(covtable) {
    out <- cbind(condition = rownames(covtable), covtable, stringsAsFactors = F)
    out <- out[!duplicated(names(out))]
    return(out)}
  covtable <- rownames_to_condition(covariates(data))

  data <- lbind(data)
  data <- base::merge(data, covtable, by = "condition", all.x = T)
  data <- as.data.frame(dplyr::filter(data, ...), stringsAsFactors = F)
  data$bloq <- ifelse(data$value <= data$lloq, "yes", "no")

  if (!is.null(transform)) data <- coordTransform(data, transform)

  if (facet == "wrap")
    p <- ggplot(data, aes(x = time, y = value, ymin = value - sigma,
                          ymax = value + sigma, group = condition, color = condition, pch = bloq)) + 
    facet_wrap(~name, scales = scales)
  if (facet == "grid")
    p <- ggplot(data, aes(x = time, y = value, ymin = value - sigma,
                          ymax = value + sigma, pch = bloq)) +  
    facet_grid(name ~ condition, scales = scales)

  p <- p + geom_point() + geom_errorbar(width = 0) +
    scale_shape_manual(name = "BLoQ", values = c(yes = 4, no = 19))

  if (all(data$bloq %in% "no"))
    p <- p + guides(shape = FALSE)

  attr(p, "data") <- data
  return(p)

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
      if (nrow(data) > 0)
        data[, C] <- condition.grid[as.character(data$condition), C]
      else
        data[, C] <- vector(mode = mode(condition.grid[[C]]))
    }
  }

  return(data)

}


#' Access the covariates in the data
#'
#' @param x Either a \link{datalist} or a \code{data.frame} with mandatory 
#' columns \code{c("name", "time", "value", "sigma", "lloq")}.
#'
#' @return The \code{condition.grid} of the data
#' @export
covariates <- function(x) {
  UseMethod("covariates", x)
}

#' @export
#' @rdname covariates
covariates.datalist <- function(x) {

  attr(x, "condition.grid")

}

#' @export
#' @rdname covariates
covariates.data.frame <- function(x) {

  exclude <- c("name", "time", "value", "sigma", "lloq")
  contains.condition <- "condition" %in% colnames(x)
  out <- unique(x[, setdiff(names(x), exclude)])

  if (contains.condition) {
    if (any(duplicated(out[["condition"]]))) {
      stop("Unique entries of condition column do not correspond to unique covariate combinations.")
    }
    rownames(out) <- out[["condition"]]
    out <- out[ , setdiff(names(out), "condition")]
  } else {
    rownames(out) <- do.call(paste_, out)
  }

  return(out)

}
