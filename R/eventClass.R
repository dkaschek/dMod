
#' @export
eventlist <- function(var = NULL, time = NULL, value = NULL, method = NULL) {
  
  out <- data.frame(var = var,
                    time = time,
                    value = value, 
                    method = method,
                    stringsAsFactors = FALSE)
  
  class(out) <- c("eventlist", "data.frame")
  return(out)
  
}

#' @export
as.eventlist <- function(x, ...) {
  UseMethod("as.eventlist", x)
}


#' @export
as.eventlist.list <- function(x) {
  
  # Check names
  required <- c("var", "time", "value", "method")
  if (!all(required %in% names(x)))
    stop("x needs to provide var, time, value and method.")
  
  # Convert to list of characters and numeric
  x <- lapply(x[required], function(element) {
    if (!is.numeric(element)) as.character(element) else element
  })
  
  # Return eventlist object
  do.call(eventlist, x)
  
  
}

#' @export
as.eventlist.data.frame <- function(x) {

  as.eventlist.list(as.list(x))
    
}

#' @export
addEvent <- function(event, ...) {
  
  UseMethod("addEvent", event)
  
}

#' @export
addEvent.eventlist <- function(event, var, time = 0, value = 0, method = "replace") {
  
  event.new <- data.frame(var = var, time = time, value = value, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}