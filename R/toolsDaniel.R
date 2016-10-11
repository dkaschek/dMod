#' Produce data frames for plotting
#' 
#' @param prediction object of class prediction list
#' @param data object of class data list
#' @param ... arguments going to subsetting of prediction and data
#' @param condition.grid data.frame with row.names according to the condition names.
#' @param substitute Substitute "_" by character \code{substitute}.
#' @return list with prediction (data.frame) and data (data.frame)
#' @export
ggdata <- function(prediction = NULL, data = NULL, ..., condition.grid = attr(data, "condition.grid"), substitute = "_") {
  
  if (!is.null(prediction)) {
    prediction <- wide2long(prediction)
    for (C in colnames(condition.grid)) {
      prediction[, C] <- condition.grid[as.character(prediction$condition), C]
    }
    prediction <- subset(prediction, ...)
    n1 <- nrow(prediction)
    
  }
    
    
  if (!is.null(data)) {
    data <- lbind(data)
    for (C in colnames(condition.grid)) {
      data[, C] <- condition.grid[as.character(data$condition), C]
    }
    data <- subset(data, ...)
  }
  
  
 
  
  out <- combine(prediction, data)
  for (n in colnames(out)) {
    mylevels <- levels(out[, n]) 
    if (!is.null(mylevels)) levels(out[, n]) <- gsub("_", substitute, mylevels, fixed = TRUE)
  }
  
  if (!is.null(prediction)) {
    return(list(prediction = out[1:n1,], data = out[-(1:n1),]))  
  } else {
    return(list(prediction = NULL, data = out))
  }
  
  
  
}

#' Produce data frames for plotting
#' 
#' @param prediction object of class prediction list
#' @param data object of class data list
#' @param ... arguments going to subsetting of prediction and data
#' @param condition.grid data.frame with row.names according to the condition names.
#' @param substitute Substitute "_" by character \code{substitute}.
#' @return list with prediction (data.frame) and data (data.frame)
#' @export
ggdata_fn <- function(prdfn = NULL, errfn = NULL, data = NULL, times, pars, ..., condition.grid = attr(data, "condition.grid"), substitute = "_") {
  
  
  
  if (!is.null(prdfn)) {
    
    prediction <- prdfn(times, pars)
    
    sigma <- NULL  
    if (!is.null(errfn)) {
      sigma <- as.prdlist(
        lapply(1:length(prediction), 
               function(i) errfn(prediction[[i]], getParameters(prediction[[i]]), conditions = names(prediction)[i])[[1]]),
        names = names(prediction)
      )
      sigma <- wide2long(sigma)
    }
    
    prediction <- wide2long(prediction)
    prediction$sigma <- NaN
    if (!is.null(sigma)) {
      common <- intersect(unique(prediction$name), unique(sigma$name))
      prediction$sigma[prediction$name %in% common] <- sigma$value[sigma$name %in% common]
    }
    
    
    for (C in colnames(condition.grid)) {
      prediction[, C] <- condition.grid[as.character(prediction$condition), C]
    }
    prediction <- subset(prediction, ...)
    n1 <- nrow(prediction)
    
  }
  
  if (!is.null(data)) {
    data <- lbind(data)
    for (C in colnames(condition.grid)) {
      data[, C] <- condition.grid[as.character(data$condition), C]
    }
    data <- subset(data, ...)
  }
  
  
  
  
  out <- combine(prediction, data)
  for (n in colnames(out)) {
    mylevels <- levels(out[, n]) 
    if (!is.null(mylevels)) levels(out[, n]) <- gsub("_", substitute, mylevels, fixed = TRUE)
  }
  
  if (!is.null(prediction)) {
    return(list(prediction = out[1:n1,], data = out[-(1:n1),]))  
  } else {
    return(list(prediction = NULL, data = out))
  }
  
  
  
}


#' @export
predict.prdfn <- function(x, times, pars, data = NULL, ...) {
  
  arglist <- list(...)
  if (any(names(arglist) == "conditions")) {
    C <- arglist[["conditions"]]
    if (!is.null(data)) {
      data <- data[C]
    }
  }
  if (is.null(data)) data <- data.frame()
  condition.grid.data <- attr(data, "condition.grid")
  
  prediction <- do.call(rbind, lapply(1:nrow(pars), function(i) {
    
    mypar <- as.parvec(pars, i)
    prediction <- x(times, mypar, deriv = FALSE, ...)
    conditions <- ifelse(is.null(names(prediction)), 1, names(prediction))
    
    condition.grid <- data.frame(row.names = conditions)
    
    # Augment by parframe metanames and obj.attributes
    mygrid <- pars[i, !colnames(pars) %in% attr(pars, "parameters")]
    mynames <- colnames(mygrid)
    if (length(mynames) > 0) {
      mynames <- paste0(".", mynames)
      colnames(mygrid) <- mynames
      condition.grid <- cbind(condition.grid, mygrid)
    }
    
    # Augment by condition.grid of data
    if (!is.null(condition.grid.data)) condition.grid <- cbind(condition.grid.data[conditions,], condition.grid)
    
    # Write condition.grid into data
    attr(data, "condition.grid") <- condition.grid

    # Return
    print(attr(data, "condition.grid"))
    as.data.frame(prediction, data = data)
    
  }))
  
  n <- nrow(prediction)
  
  if (length(data) > 0) {
    attr(data, "condition.grid") <- condition.grid.data
    data <- as.data.frame(data)
    tmp <- combine(prediction, data)
    data <- tmp[-(1:n),]
  }
  
  attr(prediction, "data") <- data
  return(prediction)  
  
  
  
}

#' Generate sample for multi-start fit
#' 
#' @param center named numeric, the center around we sample
#' @param samplefun character, indicating the random number generator,
#' defaults to \code{"rnorm"}.
#' @param fits length of the sample
#' @param ... arguments going to \code{samplefun}
#' @return matrix with the parameter samples
#' @export
mssample <- function(center, samplefun = "rnorm", fits = 20, ...) {
  
  sample.matrix <- do.call(rbind, lapply(1:fits, function(i) 
    do.call(samplefun, c(list(n = length(center), mean = center), list(...)))))
  
  colnames(sample.matrix) <- names(center)
  
  myframe <- parframe(
    cbind(
      index = 1:nrow(sample.matrix),
      as.data.frame(sample.matrix)
    ),
    metanames = "index",
    parameters= colnames(sample.matrix)
  )
  
  return(myframe)

}

