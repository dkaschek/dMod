
# Remarks --------------
# Even though this file is called dModframe"class", dMod.frame doesn't have its own class.
# The class of dMod.frames are tbl_df (read ?tibble if you want to find out more)
# I chose to call the file like this nevertheless, because I put the dMod.frame-methods, such as plotting methods in this file.



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
#' @param ... other columns, have to be named. They will be placed in a list.
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

  # To do: enlist the content of ...

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
#' @param pars A named vector of parameters to run e.g. test simulations of the model. Defaults to random parameters
#' @param times A vector of times to run e.g. test simulations of the model. Defaults to \code{seq(0, 1*t_max(data), length.out = 200)}
#'
#' @importFrom rlang quos enquo UQS
#'
#' @return The dMod.frame augmented by standardized columns
#' @export
appendObj <- function(dMod.frame,
                      prd = list(g*(x*p)),
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
#' @param parframes Expression to turn a column containing a parlist (e.g. fits) into a column of parframes
#' @param keepFits 
#' @param ... Other columns you want to mutate
#' @param keepCalls Store a record of the calls in a new colun? See \link{mutatedMod.frame}.
#'
#' @importFrom rlang quos enquo UQS
#'
#' @return The dMod.frame containing the column "parframes"
#' @export
appendParframes <- function(dMod.frame,
                            parframes = list(as.parframe(fits)),
                            keepFits = F,
                            ...,
                            keepCalls = F) {

  args <- c(list(parframes = enquo(parframes)), quos(...))

  out <- mutatedMod.frame(dMod.frame, UQS(args), keepCalls = keepCalls)
  if(!keepFits){
    out <- rowwise(out[!(names(out)=="fits")])}
  return(out)
}



# Plotting ---------------------------------

#' @export
#' @rdname plotCombined
#' @param plotErrorBands If the dMod.frame contains an observation function for the 
#' error model, use it to show error bands.
plotCombined.tbl_df <- function(model, hypothesis = 1, index = 1, ... , plotErrorBands = F) {

  dots <- substitute(alist(...))
  message("If you want to subset() the plot, specify hypothesis AND index")
  if(is.character(hypothesis)) hypothesis <- which(model$hypothesis == hypothesis)

  
  data <- model[["data"]][[hypothesis]]
  prediction <- title <- NULL
  if (is.null(model[["parframes"]])) {
    prediction <- model[["prd"]][[hypothesis]](model[["times"]][[hypothesis]], 
                                                    model[["pars"]][[hypothesis]], 
                                                    deriv = F,
                                                    fixed = model[["fixed"]][[hypothesis]])
    title <- paste(model[["hypothesis"]][[hypothesis]], "initiated with predefined (probably random) parameters")
  }
  
  if (!is.null(model[["parframes"]])) {
    myparvec <- as.parvec(model[["parframes"]][[hypothesis]], index = index)
    prediction <- model[["prd"]][[hypothesis]](model[["times"]][[hypothesis]],
                                                    pars = myparvec,
                                                    deriv = F,
                                                    fixed = model[["fixed"]][[hypothesis]])
    myvalue <- model[["parframes"]][[hypothesis]][index, "value"]
    title <- paste0(model[["hypothesis"]][[hypothesis]], "\n",
                    "value = ", round(model[["parframes"]][[hypothesis]][index,"value", drop = T],1), "\n",
                    paste0(paste(names(dots), "=", dots )[-1], collapse = "\n"))
  }


  
  myplot <- plotCombined.prdlist(prediction, data, ...) + 
    ggtitle(label = title)
  
  if (plotErrorBands) {
    predicition_with_error <- NULL
    if (!is.null(model[["e"]][[hypothesis]])){
      predicition_with_error <- dMod:::as.data.frame.prdlist(prediction, data = data, errfn = model[["e"]][[hypothesis]])
      predicition_with_error <- subset(predicition_with_error, ...)
    }
    myplot <- myplot +
      geom_ribbon(data = predicition_with_error, alpha = 0.15, size = 0)
  }  
  return(myplot)
}

#' @export
#' @rdname plotPrediction
plotPrediction.tbl_df <- function(dMod.frame, hypothesis = 1, index = 1, ... ) {

  dots <- substitute(alist(...))
  message("If you want to subset() the plot, specify hypothesis AND index")

  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)

  if (is.null(dMod.frame[["parframes"]]))
    return(
      plotPrediction.prdlist(dMod.frame[["prd"]][[hypothesis]](dMod.frame[["times"]][[hypothesis]],
                                                               dMod.frame[["pars"]][[hypothesis]],
                                                               deriv = F,
                                                               fixed = dMod.frame[["fixed"]][[hypothesis]]), ...) +
        ggtitle(paste(dMod.frame[["hypothesis"]][[hypothesis]], "initiated with predefined (probably random) parameters"))
    )


  myparvec <- as.parvec(dMod.frame[["parframes"]][[hypothesis]], index = index)
  myprediction <- dMod.frame[["prd"]][[hypothesis]](dMod.frame[["times"]][[hypothesis]],
                                                    pars = myparvec,
                                                    deriv = F,
                                                    fixed = dMod.frame[["fixed"]][[hypothesis]])
  myvalue <- dMod.frame[["parframes"]][[hypothesis]][index, "value"]

  plotPrediction.prdlist(myprediction, ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           "value = ", round(dMod.frame[["parframes"]][[hypothesis]][index,"value", drop = T],1), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}



#' @export
#' @rdname plotData
plotData.tbl_df <- function(dMod.frame, hypothesis = 1, ... ) {

  dots <- substitute(alist(...))
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)

  plotData.datalist(dMod.frame[["data"]][[hypothesis]], ...) +
    ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}



#' @export
#' @rdname plotPars
#' @param nsteps number of steps from the waterfall plot
plotPars.tbl_df <- function(dMod.frame, hypothesis = 1,  ..., nsteps = 3, tol = 1 ) {

  if (!missing(...)) {dots <- substitute(...)
  } else {
    dots <- T
  }

  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)

  if(length(hypothesis)==1) {
    myparframe <- dMod.frame[["parframes"]][[hypothesis]]
    myparframe <- getSteps(myparframe, nsteps = nsteps, tol = tol)

    plotPars.parframe(myparframe, ..., tol = tol) +
      ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                             "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                             paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
  } else {
    mydata <- do.call(rbind,lapply(hypothesis, function(hyp) {
      pars <- dMod.frame[["parframes"]][[hyp]]
      pars <- cbind(index = 1:nrow(pars), pars[order(pars$value),!(names(pars) == "index")])
      steps <- getStepIndices(pars, nsteps, tol)
      pars <- getSteps(pars, nsteps = nsteps, tol = tol)
      cbind(hypothesis = hyp,step = factor(findInterval(pars$index, steps), labels = steps), pars)
    }))
    mydata$hypothesis <- as.factor(mydata$hypothesis)

    myindices <- with(mydata, eval(dots))
    mydata <- mydata[myindices,]

    mydata2 <- wide2long.data.frame(mydata[-c(4,5,6)], keep = 1:3)
    mydata2$hyp_step <- paste(mydata2$hypothesis, mydata2$step, sep = "_")

    ggplot2::ggplot(mydata2, aes(x = name, y = value, color = hyp_step)) +
      geom_boxplot(outlier.alpha = 0) +     # annotate("text", x = jumps + 1, y = y.jumps, label = jumps, hjust = 0, color = "red", size = 3) +
      xlab("parameter") + ylab("value") + theme_dMod() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  }

}

#' @export
#' @rdname plotValues
plotValues.tbl_df <- function(dMod.frame, hypothesis = 1, ..., tol = 1 ) {

  dots <- substitute(...)

  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)

  if(length(hypothesis)==1) {
    plotValues.parframe(dMod.frame[["parframes"]][[hypothesis]], tol = tol, ...) +
      ggtitle(label = paste0(dMod.frame[["hypothesis"]][[hypothesis]], "\n",
                             "best value = ", round(dMod.frame[["parframes"]][[hypothesis]][1,"value", drop = T],1), "\n",
                             paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
  } else {
    mydata <- do.call(rbind,lapply(hypothesis, function(hyp) {
      pars <- dMod.frame[["parframes"]][[hyp]]
      pars <- cbind(index = 1:nrow(pars), pars[order(pars$value),])

      cbind(hypothesis = hyp, pars)
    }))
    mydata$hypothesis <- as.factor(mydata$hypothesis)

    myindices <- with(mydata, eval(dots))
    mydata <- mydata[myindices,]

    ggplot2::ggplot(mydata, aes(x = index, y = value, pch = converged, color = hypothesis)) +
      # geom_vline(xintercept = jumps, lty = 2) +
      geom_point() +
      # annotate("text", x = jumps + 1, y = y.jumps, label = jumps, hjust = 0, color = "red", size = 3) +
      xlab("index") + ylab("value") + theme_dMod()

  }
}



#' @export
#' @rdname plotProfile
#' @param hypothesis numeric, can be longer than 1
plotProfile.tbl_df <- function(dMod.frame, hypothesis = 1, ...) {
  dots <- substitute(alist(...))

  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)

  myprofs <- dMod.frame[["profiles"]][hypothesis] %>% setNames(dMod.frame[["hypothesis"]][hypothesis])

  plotProfile.list(myprofs, ...) +
    ggtitle(label = paste0("hypotheses: ", paste0(hypothesis, collapse = ", "), "\n",
                           "best values = ", paste0(lapply(hypothesis, function(hyp) round(dMod.frame[["parframes"]][[hyp]][1,"value", drop = T],1)), collapse = ", "), "\n",
                           paste0(paste(names(dots), "=", dots )[-1], collapse = "\n")) )
}

#' @export
#' @rdname plotPaths
#' @param hypothesis numeric, can be longer than 1
plotPaths.tbl_df <- function(dMod.frame, hypothesis = 1, ...) {
  dots <- substitute(alist(...))
  
  if(is.character(hypothesis)) hypothesis <- which(dMod.frame$hypothesis == hypothesis)
  
  myprofs <- dMod.frame[["profiles"]][hypothesis] %>% setNames(dMod.frame[["hypothesis"]][hypothesis])
  
  plotPaths.list(myprofs, ...)
}


# Other methods -----

#' @export 
#' @param hypothesis The hypothesis in the dMod.frame
#' @rdname covariates
covariates.tbl_df <- function(x, hypothesis = 1) {
  covariates.datalist(x[["data"]][[hypothesis]])
}



# Get functions --------------

#' @export
#' @param hypothesis The hypothesis in the dMod.frame
#' @rdname getObservables
getObservables.tbl_df <- function(x, hypothesis = 1, ...) {
  dMod.frame <- x
  Reduce(union, lapply(getEquations(dMod.frame[["g"]][[hypothesis]]), names))
}


#' @export
#' @param hypothesis The hypothesis in the dMod.frame
#' @rdname getConditions
getConditions.tbl_df <- function(x, hypothesis = 1, ...) {
  dMod.frame <- x
  getConditions.fn(dMod.frame$obj[[hypothesis]])
}



# Interaction with .GlobalEnv ----

#' Load one row of a dMod.frame into the .GlobalEnv
#'
#' @param dMod.frame A dMod.frame
#' @param hypothesis character or numeric. specifying the name  or the index of the hypothesis
#' @param prefix Prefix appended to the object names in .GlobalEnv
#' @param suffix Suffix appended to the object names in .GlobalEnv
#'
#' @export
#'
#' @examples
#' testframe <- dplyr::tibble(hypothesis = c("linear", "quadratic"),
#'                     plots = list(plot(1:10,1:10), plot(1:10,(1:10)^2)),
#'                     myfun = list(function(x,a) {a * x}, function(x,a) {a * x^2}),
#'                     a = c(1:2))
#'
#'
#' checkout_hypothesis(testframe, "quadratic", prefix = "quad")
#' quadplots
#' quadmyfun(1:10, quada)
checkout_hypothesis <- function(dMod.frame, hypothesis, prefix = "", suffix = "") {

  if(is.numeric(hypothesis)) {
    mydMod.frame <- dMod.frame[hypothesis,]
  } else {
    mydMod.frame <- dMod.frame[dMod.frame[["hypothesis"]]==hypothesis,]
  }
  lapply(seq_along(mydMod.frame), function(i)  {
    value <- mydMod.frame[[i]]
    if(is.list(value)&length(value)==1) value = value[[1]]
    try(assign(x = paste0(prefix,names(mydMod.frame)[i],suffix),
               value = value,
               pos = .GlobalEnv),
        silent = T)
  })

  message("Objects in .GlobalEnv are updated")

}












