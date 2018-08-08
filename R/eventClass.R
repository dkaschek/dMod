#' Eventlist
#' 
#' An eventlist is a data.frame with the necessary parameters to define an event as columns and specific events as rows.
#' Event time and value can be passed as parameters, which can also be estimated.
#' 
#' The function \code{addEvent} is pipe-friendly
#'  
#' @param var Character, the state to which the event is applied
#' @param time Character or Numeric, the time at which the event happens
#' @param value Character or Numeric, the value of the event
#' @param method Character, options are "replace", "add" or "multiply"
#'
#' @return data.frame with class eventlist
#' @export
#'
#' @examples
#' eventlist(var = "A", time = "5", value = 1, method = "add")
#' 
#' events <- addEvent(NULL, var = "A", time = "5", value = 1, method = "add")
#' events <- addEvent(events, var = "A", time = "10", value = 1, method = "add")
eventlist <- function(var = NULL, time = NULL, value = NULL, method = NULL) {
  
  out <- data.frame(var = var,
                    time = time,
                    value = value, 
                    method = method,
                    stringsAsFactors = FALSE)
  
  class(out) <- c("eventlist", "data.frame")
  return(out)
  
}

#' Coerce to eventlist
#'
#' @param ... not used
#' @export
as.eventlist <- function(x, ...) {
  UseMethod("as.eventlist", x)
}


#' @export
#' @rdname as.eventlist
#' @param x list, data.frame
as.eventlist.list <- function(x, ...) {
  
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
#' @rdname as.eventlist
as.eventlist.data.frame <- function(x, ...) {

  as.eventlist.list(as.list(x))
    
}

#' @rdname eventlist
#' @param event object of class \code{eventlist}
#' @param ... not used
#' @export
addEvent <- function(event, ...) {
  
  UseMethod("addEvent", event)
  
}

#' @export
addEvent.eventlist <- function(event, var, time = 0, value = 0, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}

#' @export
addEvent.NULL <- function(event, var, time = 0, value = 0, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}