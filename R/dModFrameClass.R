# dMod.frame constructor ----------------------------------------------------------------


#' Generate a dMod.frame
#'
#' @description Basically, a dMod.frame is a \link{tibble}, which is grouped \link{rowwise}.
#' 
#' Since the dMod.frame is also designed for interactive use, the class will not be called
#' "dMod.frame" as I initially planned, but will be c("tbl_df").
#' This way, I don't have to overwrite all the dplyr-verbs.
#' 
#' The dMod.frame object stores all objects that are needed to reproducibly
#' produce a result, such as a plot or profiles, in one row of a \link{tibble}.
#' Each row corresponds to a distinct hypothesis, e.g. one would have two distinct rows
#' for fitting with and without a prior.
#'
#'
#'
#' @param hypothesis Character. Description of the hypothesis
#' @param g fn
#' @param x fn
#' @param p fn
#' @param data data.frame or datalist, will be coerced to datalist
#' @param e fn
#' @param ...
#'
#' @return Object of class \code{tbl_df}, which is grouped rowwise.
#'
#' @importFrom dplyr tibble rowwise
#'
#' @export
#'
#' @example inst/examples/dMod.frame.R
dMod.frame <- function(hypothesis, g, x, p, data, e = NULL,...) {

  enlist <- function(x) {
    if (is.list(x) & (!(is.data.frame(x)|is.datalist(x))) ) return(x)
    else return(list(x))
  }

  out <- tibble(hypothesis = hypothesis,
         g = enlist(g),
         x = enlist(x),
         p = enlist(p),
         data = enlist(as.datalist(data)),
         e = enlist(e),
         ...)
  out <- rowwise(out)

  return(out)
}


#' A version of dplyr::mutate
#'
#' @description This is basically dplyr::mutate with two differences
#' 1. It groups the tibble \link{rowwise} before mutating
#' 2. It allows to store the calls. This is intended for use when your objects are
#' not standard transformations, such as \code{prd = g*x*p} but more complicated and 
#' you want to keep a record of what you did.
#' 
#' If the result of your ... expressions is not atomic, make sure to wrap your 
#' expression in \code{list()}
#'
#' @param dMod.frame A dMod.frame
#' @param ... Expressions going to mutate()
#' @param keepCalls Should the dots ... be recorded?
#'
#' @importFrom dplyr rowwise mutate
#' @importFrom rlang quos UQS
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mytbl <- tibble(a = 1:2, b = letters[1:2]) %>%
#' mutatedMod.frame(e = paste0(a,b), keepCalls = T) %>%
#' mutatedMod.frame(f = paste0(e, "=", a, "*", b), keepCalls = T) %>%
#' mutatedMod.frame(e = paste0(a,"*", b), keepCalls = T)
#'
#' mytbl$.calls
#' }
mutatedMod.frame <- function(dMod.frame,
                             ...,
                             keepCalls = F) {

  args <- quos(...)

  if (keepCalls) {
    if (is.null(dMod.frame[[".calls"]])) dMod.frame <- mutate(dMod.frame, .calls = list(NULL))
    return(mutate(rowwise(dMod.frame), UQS(args), .calls = list(c(.calls, list(args)))))
  } else {
    return(mutate(rowwise(dMod.frame), UQS(args)))
  }

}


#' Append an objective function to a basic dMod.frame
#'
#' @param dMod.frame A dMod.frame
#' @param prd Expression after which the concatenated prediction function is formed. Has to wrapped in list()
#' @param obj_data Expression after which the objective function which describes the data is formed. Has to wrapped in list()
#' @param obj This object is taken by the standard fitting functions. At typical expression would be \code{list(obj_data +  constr)}. Has to wrapped in list()
#' @param ... Other columns which are mutations of existing ones or new columns.
#' @param keepCalls Store a record of the calls in a new colun? See \link{mutatedMod.frame}.
#'
#' @importFrom rlang quos enquo UQS
#'
#' @return The dMod.frame augmented by standardized columns
#' @export
appendObj <- function(dMod.frame,
                      prd = list(g*x*p),
                      obj_data = list(normL2(data, prd, e)),
                      obj = list(obj_data),
                      pars = list(structure(rnorm(length(getParameters(obj))), names = getParameters(obj))),
                      times = list(seq(min(as.data.frame(data)[["time"]]), max(as.data.frame(data)[["time"]])*1.1, length.out = 200)),
                      ...,
                      keepCalls = F) {

  args <- c(list(prd = enquo(prd),
                 obj_data = enquo(obj_data),
                 obj = enquo(obj),
                 pars = enquo(pars),
                 times = enquo(times)
                 ),
            quos(...))

  mutatedMod.frame(dMod.frame, UQS(args), keepCalls = keepCalls)

}



#' Make a column "parframes" out of "fits"
#'
#' Most plotting functions rely on a column "parframes" to be existent in the dMod.frame
#'
#' @param dMod.frame A dmod.frame, preferably with a column \code{fits}.
#' @param parframes Expression to turn a column containing a parlist (e.g. fits) into a column of parframess
#' @param ... Other columns you want to mutate.
#' @param keepCalls Store a record of the calls in a new colun? See \link{mutatedMod.frame}.
#' 
#' @importFrom rlang quos enquo UQS
#'
#' @return The dMod.frame containing the column "parframes"
#' @export
appendParframes <- function(dMod.frame,
                            parframes = list(as.parframe(fits)),
                            ...,
                            keepCalls = F) {

  args <- c(list(parframes = enquo(parframes)), quos(...))

  mutatedMod.frame(dMod.frame, UQS(args), keepCalls = keepCalls)

}



# Plotting ---------------------------------

#' @rdname plotCombined.prdlist
#' 
#' PlotCombined method for dMod.frame
#'
#' @param dMod.frame A dMod.frame
#' @param hypothesis index specifying the row (hypothesis)
#' @param index index specifying the index of the fit going to \code{parframes %>% as.parvec(index)}
#' @param ... Arguments going to subset
#'
#' @export
plotCombined.tbl_df <- function(dMod.frame, hypothesis = 1, index = 1, ... ) {
  
  dots <- substitute(alist(...))
  
  message("If you want to subset() the plot, specify hypothesis and index")
  
  
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  # i <- hypothesis #so i can copy other code
  
  times <- NULL
  if (!is.null(dMod.frame[["times"]]))
    times <- dMod.frame[["times"]][[hypothesis]]
  else {
    times <- as.data.frame(dMod.frame[["data"]][[hypothesis]])[["time"]]
    times <- seq(min(times), max(times)*1.1, length.out = 100)
  }
  
  if (is.null(dMod.frame[["parframes"]]))
    return(
      plotCombined.prdlist(
        dMod.frame[["prd"]][[hypothesis]](times, dMod.frame[["pars"]][[hypothesis]], deriv = F),
        dMod.frame[["data"]][[hypothesis]],
        ...) +
        ggtitle(paste(dMod.frame[["hypothesis"]][[hypothesis]], "initiated with predefined (probably random) parameters"))
    )
  
  
  
  myparvec <- as.parvec(dMod.frame[["parframes"]][[hypothesis]], index = index)
  
  myprediction <- dMod.frame[["prd"]][[hypothesis]](times,
                                                    pars = myparvec,
                                                    deriv = F)
  
  myvalue <- dMod.frame[["parframes"]][[hypothesis]][index, "value"]
  
  plotCombined.prdlist(myprediction,  dMod.frame[["data"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "value = ", round(dMod.frame[["parframes"]][[hypothesis]][index,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}





#' @export
#' @rdname plotData
plotData.tbl_df <- function(dMod.frame, hypothesis = 1, ... ) {
  
  dots <- substitute(alist(...))
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  # i <- hypothesis #so i can copy other code
  
  # myparvec <- as.parvec(dMod.frame[["parframes"]][[hypothesis]], index = index)
  
  # myprediction <- dMod.frame[["prd"]][[hypothesis]](times = seq(0, max(as.data.frame(dMod.frame[["data"]][[hypothesis]])$time),length.out = 100),
  #                                                   pars = myparvec,
  #                                                   deriv = F)
  
  # myvalue <- dMod.frame[["parframes"]][[hypothesis]][1, "value"]
  
  plotData.datalist(dMod.frame[["data"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}



#' @export
#' @rdname plotPars
plotPars.tbl_df <- function(dMod.frame, hypothesis = 1, ... ) {
  
  dots <- substitute(alist(...))
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  
  plotPars.parframe(dMod.frame[["parframes"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}

#' @export
#' @rdname plotValues
plotValues.tbl_df <- function(dMod.frame, hypothesis = 1, ... ) {
  
  dots <- substitute(alist(...))
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  
  plotValues.parframe(dMod.frame[["parframes"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}



#' @export
#' @rdname plotProfile
plotProfile.tbl_df <- function(dMod.frame, hypothesis = 1, ...) {
  dots <- substitute(alist(...))
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  
  plotProfile.list(dMod.frame[["profiles"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}
