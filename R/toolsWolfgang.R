#' Pad string to desired width
#'
#' @param string String to pad
#' @param width Desired width of padded string
#' @param where Padding can be inserted to the right or left of <string>.
#'        Default to 'right'.
#' @param padding A single character with with the padding space is filled.
#'        Defaults to blank ' ' yielding invisible padding.
#' @param autoelide If TRUE, <string> is elided if it is wider than <width>. The
#'        position of eliding follows <where>. Defaults to FALSE.
#'
#' @return Padded string of length <width>.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
strpad <- function(string, width, where = "right", padding = " ", autoelide = FALSE) {

  # Check function arguments
  # <padding> is assumed to have exactry one character
  if (nchar(padding) != 1) {
    stop("<padding> must be of width 1")
  }

  # We only work on a character vector
  if (!inherits(string, "character")) {
    warning("Argument <string> does not inherit class 'character'")
    return(NULL)
  }


  # Assemble parameters
  strWidth <- nchar(string)


  # Check if eliding is needed
  if (strWidth > width) {
    if (!autoelide) {
      stop("<string> too wide. Consider autoelide = TRUE")

    } else {
      return(strelide(string, width, where))
    }
  }


  # Do the padding
  m_padding <- paste0(rep(padding, width - strWidth), collapse = "")
  if (where == "left") {
    return(paste0(m_padding, string))

  } else if (where == "right") {
    return(paste0(string, m_padding))

  } else {
    stop("<where> must be 'left' or 'right'")
  }
}



#' Elide character vector
#'
#' @param string String subject to eliding
#' @param width Width including eliding ... of return string
#' @param where Eliding can happen at 'left', 'middel', or 'right'. Defaults to
#'        'right'.
#' @param force Elide, even is <string> is shorter than <width>. Default to
#'        'FALSE'.
#'
#' @return Elided string of length <width>.
#'
#' @details
#' Elide a string to <width>. Eliding can happen at 'left', 'middle', or
#' 'right'. #' If forcing = FALSE, which is the default, strings shorten than
#' <width> are returend unaltered; forcing = TRUE inserts eliding symbols (...)
#' in any case.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
strelide <- function(string, width, where = "right", force = FALSE) {

  # Functions for eliding
    # Eliding in the middle
  elideMiddle <- function(string, width) {
    if (width == 1) {
      return('.')
    }
    strWidth <- nchar(string)
    frontWidth <- max(ceiling((width - 3) / 2), 1)
    endWitdh <- floor((width - 3) / 2)
    string <- paste0(substr(string, 1, frontWidth),  "...",
                     substr(string, strWidth - endWitdh + 1, strWidth))
    return(strtrim(string, width))
  }

  # Elide on the right
  elideRight <- function(string, width) {
    if (width <= 3) {
      return(elideMiddle(string, width))
    } else {

      substr(string, width - 2, width) <- "..."
      return(strtrim(string, width))
    }
  }

  # Eliding on the left
  elideLeft <- function(string, width) {
    string <- paste(rev(substring(string, 1:nchar(string), 1:nchar(string))), collapse="")
    string <- elideRight(string, width)
    return(paste(rev(substring(string, 1:nchar(string), 1:nchar(string))), collapse=""))
  }

  # Check function argumets
  # We only work on a character vector
  if (!inherits(string, "character")) {
    warning("Argument does not inherit class 'character'")
    return(NULL)
  }


  # Assemble parameters
  strWidth <- nchar(string)


  # Eliding
  # Do we need to elide?
  if (strWidth <= width && !force) {
    return(string)
  }


  # Are we forced to elide?
  if (width > strWidth && force) {
    width <- strWidth
  }


  # Eliding the normal case
  if (where == "left") {
    return(elideLeft(string, width))

  } else if (where == "middle") {
    return(elideMiddle(string, width))

  } else if (where == "right") {
    return(elideRight(string, width))

  } else {
    stop("<where> must be 'left', 'middle', or 'right'")
  }
}




#' Non-Linear Optimization, multi start
#' 
#' @description Wrapper around \code{\link{trust}} allowing for multiple fits 
#'   from randomly chosen initial values.
#'   
#' @param objfun Objective function, see \code{\link{trust}}.
#' @param center Parameter values around which the initial values for each fit 
#'   are randomly sampled. The initial values handed to \link{trust} are the sum
#'   of center and the output of \option{samplefun}, center + 
#'   \option{samplefun}. See \code{\link{trust}}, parinit.
#' @param studyname The names of the study or fit. This name is used to 
#'   determine filenames for interim and final results. See Details.
#' @param rinit Starting trust region radius, see \code{\link{trust}}.
#' @param rmax Maximum allowed trust region radius, see \code{\link{trust}}.
#' @param fits Number of fits (jobs).
#' @param cores Number of cores for job parallelization.
#' @param samplefun Function to sample random initial values. It is assumed, 
#'   that \option{samplefun} has a named parameter "n" which defines how many 
#'   random numbers are to be returned, such as for \code{\link{rnorm}} or 
#'   \code{\link{runif}}. By default \code{\link{rnorm}} is used. Parameteres 
#'   for samplefun are simply appended as named parameters to the mstrust call 
#'   and automatically handed to samplefun by matching parameter names.
#' @param resultPath character indicating the folder where the results should 
#'   be stored. Defaults to ".". 
#' @param stats If true, the same summary statistic as written to the logfile is
#'   printed to command line on mstrust completion.
#' @param ... Additional parameters handed to trust(), samplefun(), or the 
#'   objective function by matching parameter names. All unmatched parameters 
#'   are handed to the objective function objfun(). The log file starts with a 
#'   table telling which parameter was assigend to which function.
#'   
#' @details By running multiple fits starting at randomly chosen inital 
#'   parameters, the chisquare landscape can be explored using a deterministic 
#'   optimizer. Here, \code{\link{trust}} is used for optimization. The standard
#'   procedure to obtain random initial values is to sample random variables 
#'   from a uniform distribution (\code{\link{rnorm}}) and adding these to 
#'   \option{center}. It is, however, possible, to employ any other sampling 
#'   strategy by handing the respective function to mstrust(), 
#'   \option{samplefun}.
#'   
#'   In case a special sampling is required, a customized sampling function can 
#'   be used. If, e.g., inital values leading to a non-physical systems are to 
#'   be discarded upfront, the objective function can be addapted accordingly.
#'   
#'   All started fits either lead to an error or complete converged or
#'   unconverged. A statistics about the return status of fits can be shown by
#'   setting \option{stats} to TRUE.
#'   
#'   Fit final and intermediat results are stored under \option{studyname}. For
#'   each run of mstrust for the same study name, a folder under
#'   \option{studyname} of the form "trial-x-date" is created. "x" is the number
#'   of the trial, date is the current time stamp. In this folder, the
#'   intermediate results are stored. These intermediate results can be loaded
#'   by \code{\link{load.parlist}}. These are removed on successfull completion
#'   of mstrust. In this case, the final list of fit parameters
#'   (parameterList.Rda) and the fit log (mstrust.log) are found instead.
#'   
#' @return A parlist holding errored and converged fits.
#'   
#' @seealso \code{\link{trust}}, \code{\link{rnorm}}, \code{\link{runif}}, 
#'   \code{\link{as.parframe}}
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'  
#' @example inst/examples/test_blocks.R
#'     
#' @export
mstrust <- function(objfun, center, studyname, rinit = .1, rmax = 10, fits = 20, cores = 1,
                    samplefun = "rnorm", resultPath = ".", stats = FALSE,
                    ...) {

  narrowing <- NULL
  
  # Argument parsing, sorting, and enhancing
  # Gather all function arguments
  varargslist <- list(...)

  argslist <- as.list(formals())
  argslist <- argslist[names(argslist) != "..."]

  argsmatch <- as.list(match.call(expand.dots = TRUE))
  namesinter <- intersect(names(argslist), names(argsmatch))

  argslist[namesinter] <- argsmatch[namesinter]
  argslist <- c(argslist, varargslist)

  # Add extra arguments
  argslist$n <- length(center) # How many inital values do we need?

  # Determine target function for each function argument.
  # First, define argument names used locally in mstrust().
  # Second, check what trust() and samplefun() accept and check for name clashes.
  # Third, whatever is unused is passed to the objective function objfun().
  nameslocal <- c("studyname", "center", "fits", "cores", "samplefun",
                  "resultPath", "stats", "narrowing")
  namestrust <- intersect(names(formals(trust)), names(argslist))
  namessample <- intersect(names(formals(samplefun)), names(argslist))
  if (length(intersect(namestrust, namessample) != 0)) {
    stop("Argument names of trust() and ", samplefun, "() clash.")
  }
  namesobj <- setdiff(names(argslist), c(namestrust, namessample, nameslocal))


  # Assemble argument lists common to all calls in mclapply
  # Sample function
  argssample <- structure(vector("list", length = length(namessample)), names = namessample)
  for (name in namessample) {
    argssample[[name]] <- argslist[[name]]
  }

  # Objective function
  argsobj <- structure(vector("list", length = length(namesobj)), names = namesobj)
  for (name in namesobj) {
    argsobj[[name]] <- argslist[[name]]
  }

  # Trust optimizer, except for initial values
  argstrust <- structure(vector("list", length = length(namestrust)), names = namestrust)
  for (name in namestrust) {
    argstrust[[name]] <- argslist[[name]]
  }


  # Assemble and create output filenames, folders and files
  m_timeStamp <- paste0(format(Sys.time(), "%d-%m-%Y-%H%M%S"))
  
  # Folders
  resultFolderBase <- file.path(argslist$resultPath, argslist$studyname)
  m_trial <- paste0("trial-", length(dir(resultFolderBase, pattern = "trial*")) + 1)
  resultFolder <- file.path(resultFolderBase, paste0(m_trial, "-", m_timeStamp))
  
  interResultFolder <- file.path(resultFolder, "interRes")
  dir.create(path = interResultFolder, showWarnings = FALSE, recursive = TRUE)
  
  # Files
  fileNameLog <- paste0("mstrust.log")
  fileNameParList <- paste0("parameterList.Rda")
  fileLog <- file.path(resultFolder, fileNameLog)
  fileParList <- file.path(resultFolder, fileNameParList)
  
  
  
  # Apply trust optimizer in parallel
  # The error checking leverages that mclappy runs each job in a try().
  logfile <- file(fileLog, open = "w")
  
  # Parameter assignment information
  if (is.null(narrowing) || narrowing[1] == 1) {
    msg <- paste0("Parameter assignment information\n",
                  strpad("mstrust", 12),                        ": ", paste0(nameslocal, collapse = ", "), "\n",
                  strpad("trust", 12),                          ": ", paste0(namestrust, collapse = ", "), "\n",
                  strpad(as.character(argslist$samplefun), 12), ": ", paste0(namessample, collapse = ", "), "\n\n")
                  #strpad(as.character(argslist$objfun), 12),    ": ", paste0(namesobj, collapse = ", "), "\n\n")
    writeLines(msg, logfile)
    flush(logfile)
  }

  # Write narrowing status information to file
  if (!is.null(narrowing)) {
    msg <- paste0("--> Narrowing, run ", narrowing[1], " of ", narrowing[2], "\n",
                  "--> " , fits, " fits to run\n")
    writeLines(msg, logfile)
    flush(logfile)
  }

  m_parlist <- as.parlist(parallel::mclapply(1:fits, function(i) {
    argstrust$parinit <- center + do.call(samplefun, argssample)
    fit <- do.call(trust, c(argstrust, argsobj))

    # Keep only numeric attributes of object returned by trust()
    attr.fit <- attributes(fit)
    keep.attr <- sapply(attr.fit, is.numeric)
    fit <- fit[1:length(fit)] # deletes attributes
    if (any(keep.attr)) attributes(fit) <- c(attributes(fit), attr.fit[keep.attr]) # attach numeric attributes
    
    
    # In some crashes a try-error object is returned which is not a list. Since
    # each element in the parlist is assumed to be a list, we wrap these cases.
    if (!is.list(fit)) {
      f <- list()
      f$error <- fit
      fit <- f
    }
    
    fit$parinit <- argstrust$parinit

    # Write current fit to disk
    saveRDS(fit, file = file.path(interResultFolder, paste0("fit-", i, ".Rda")))

    # Reporting
    # With concurent jobs and everyone reporting, this is a classic race
    # condition. Assembling the message beforhand lowers the risk of interleaved
    # output to the log.
    msgSep <- "-------"
    if (any(names(fit) == "error")) {
      msg <- paste0(msgSep, "\n",
                    "Fit ", i, " failed after ", fit$iterations, " iterations with error\n",
                    "--> ", fit$error,
                    msgSep, "\n")

      writeLines(msg, logfile)
      flush(logfile)
    } else {
      msg <- paste0(msgSep, "\n",
                    "Fit ", i, " completed\n",
                    "--> iterations : ", fit$iterations, "\n",
                    "-->  converged : ", fit$converged, "\n",
                    "--> obj. value : ", round(fit$value, digits = 2), "\n",
                    msgSep)

      writeLines(msg, logfile)
      flush(logfile)
    }

    return(fit)
  }, mc.preschedule = FALSE, mc.silent = FALSE, mc.cores = cores))
  close(logfile)


  # Cull failed and completed fits Two kinds of errors occure. The first returns
  # an object of class "try-error". The reason for these failures are unknown to
  # me. The second returns a list of results from trust(), where one name of the
  # list is error holding an object of class "try-error". These abortions are 
  # due to errors which are captured within trust(). Completed fits return with 
  # a valid result list from trust(), with "error" not part of its names. These
  # fits, can still be unconverged, if the maximim number of iterations was the
  # reason for the return of trust(). Be also aware of fits which converge due
  # to the trust radius hitting rmin. Such fits are reported as converged but
  # are not in truth.
  m_trustFlags.converged = 0
  m_trustFlags.unconverged = 1
  m_trustFlags.error = 2
  m_trustFlags.fatal = 3
  idxStatus <- sapply(m_parlist, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error")) {
      return(m_trustFlags.error)
    } else if (!any(names(fit) == "converged")) {
      return(m_trustFlags.fatal)
    } else if (fit$converged) {
      return(m_trustFlags.converged)
    } else {
      return(m_trustFlags.unconverged)
    }
  })

  
  # Wrap up
  # Write out results
  saveRDS(m_parlist, file = fileParList)

  # Remove temporary files
  unlink(interResultFolder, recursive = TRUE)
  

  # Show summary
  sum.error <- sum(idxStatus == m_trustFlags.error)
  sum.fatal <- sum(idxStatus == m_trustFlags.fatal)
  sum.unconverged <- sum(idxStatus == m_trustFlags.unconverged)
  sum.converged <- sum(idxStatus == m_trustFlags.converged)
  msg <- paste0("Mutli start trust summary\n",
                "Outcome     : Occurrence\n",
                "Error       : ", sum.error, "\n",
                "Fatal       : ", sum.fatal, " must be 0\n",
                "Unconverged : ", sum.unconverged, "\n",
                "Converged   : ", sum.converged, "\n",
                "           -----------\n",
                "Total       : ", sum.error + sum.fatal + sum.unconverged + sum.converged, paste0("[", fits, "]"), "\n")
  logfile <- file(fileLog, open = "a")
  writeLines(msg, logfile)
  flush(logfile)
  close(logfile)

  if (stats) {
    cat(msg)
  }

  return(m_parlist)
}




#' Construct fitlist from temporary files.
#'
#' @description An aborted \code{\link{mstrust}}
#'   leaves behind results of already completed fits. This command loads these
#'   fits into a fitlist.
#'
#' @param folder Path to the folder where the fit has left its results.
#'
#' @details The commands \code{\link{mstrust}} save
#'   each completed fit along the multi-start sequence such that the results can
#'   be resurected on abortion. This command loads a fitlist from these
#'   intermediate results.
#'
#' @return An object of class parlist.
#'
#' @seealso \code{\link{mstrust}}
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
load.parlist <- function(folder) {
  # Read in all fits
  m_fileList <- dir(folder, pattern = "*.Rda")
  m_parVec <- lapply(m_fileList, function(file) {
    return(readRDS(file.path(folder, file)))
  })
  
  return(as.parlist(m_parVec))
}






#' Select a parameter vector from a parameter frame.
#' 
#' @description Obtain a parameter vector from a parameter frame.
#' 
#' @param x A parameter frame, e.g., the output of
#'   \code{\link{as.parframe}}.
#' @param index Integer, the parameter vector with the \code{index}-th lowest
#'   objective value.
#' @param ... not used right now
#'   
#' @details With this command, additional information included in the parameter
#'   frame as the objective value and the convergence state are removed and a
#'   parameter vector is returned. This parameter vector can be used to e.g.,
#'   evaluate an objective function.
#'   
#'   On selection, the parameters in the parameter frame are ordered such, that
#'   the parameter vector with the lowest objective value is at \option{index}
#'   1. Thus, the parameter vector with the \option{index}-th lowest objective
#'   value is easily obtained.
#'   
#' @return The parameter vector with the \option{index}-th lowest objective
#'   value.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
as.parvec.parframe <- function(x, index = 1, ...) {
  parframe <- x
  m_order <- order(parframe$value)
  metanames <- attr(parframe, "metanames")
  best <- as.parvec(unlist(parframe[m_order[index], attr(parframe, "parameters")]))
  if ("converged" %in% metanames && !parframe[m_order[index],]$converged) {
    warning("Parameter vector of an unconverged fit is selected.", call. = FALSE)
    }
  return(best)
}



#' Select attributes.
#' 
#' @description Select or discard attributes from an object.
#'   
#' @param x The object to work on
#' @param atr An optional list of attributes which are either kept or removed. 
#'   This parameter defaults to dim, dimnames, names,  col.names, and row.names.
#' @param keep For keep = TRUE, atr is a positive list on attributes which are 
#'   kept, for keep = FALSE, \option{atr} are removed.
#'   
#' @return x with selected attributes.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
attrs <- function(x, atr = NULL, keep = TRUE) {

  if (is.null(atr)) {
    atr <- c("class", "dim", "dimnames", "names", "col.names", "row.names")
  }
  
  xattr <- names(attributes(x))
  if (keep == TRUE) {
    attributes(x)[!xattr %in% atr] <- NULL
  } else {
    attributes(x)[xattr %in% atr] <- NULL
  }
  
  return(x)
}



#' Print object and its "default" attributes only.
#' 
#' @param x Object to be printed
#' @param list_attributes Prints the names of all attribute of x, defaults to 
#'   TRUE
#'   
#' @details Before the \option{x} is printed by print.default, all its arguments
#'   not in the default list of \code{\link{attrs}} are removed.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, 
#'   \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
print0 <- function(x, list_attributes = TRUE ) {
  if (list_attributes == TRUE) {
    cat("List of all attributes: ", names(attributes(x)), "\n")
  }
  
  print.default(attrs(x))
}



#' Reduce replicated measurements to mean and standard deviation
#'
#' @description
#' Obtain the mean and standard deviation from replicates per condition.
#'
#' @param file Data file of csv. See Format for details.
#' @param select Names of the columns in the data file used to define
#'        conditions, see Details.
#' @param datatrans Character vector describing a function to transform data.
#'        Use \kbd{x} to refere to data.
#'
#'
#' @format
#' The following columns are mandatory for the data file.
#' \describe{
#'  \item{name}{Name of the observed species.}
#'  \item{time}{Measurement time point.}
#'  \item{value}{Measurement value.}
#'  \item{condition}{The condition under which the observation was made.}
#' }
#'
#' In addition to these columns, any number of columns can follow to allow a
#' fine grained definition of conditions. The values of all columns named in
#' \option{select} are then merged to get the set of conditions.
#'
#' @details
#' Experiments are usually repeated multiple times possibly under different
#' conditions leading to replicted measurements. The column "Condition" in the
#' data allows to group the data by their condition. However, sometimes, a more
#' fine grained grouping is desirable. In this case, any number of additional
#' columns can be append to the data. These columns are referred to as
#' "condition identifier". Which of the condition identifiers are used to do the
#' grouping is user defined by anouncing the to \option{select}. The mandatory
#' column "Condition" is always used. The total set of different conditions is
#' thus defined by all combinations of values occuring in the selected condition
#' identifiers. The replicates of each condition are then reduced to mean and
#' variance.New conditions names are derived by merging all conditions which
#' were used in mean and std.
#'
#' @return
#' A data frame of the following variables
#' \describe{
#'  \item{time}{Measurement time point.}
#'  \item{name}{Name of the observed species.}
#'  \item{value}{Mean of replicates.}
#'  \item{sigma}{Standard error of the mean, NA for single measurements.}
#'  \item{n}{The number of replicates reduced.}
#'  \item{condition}{The condition for which the value and sigma were calculated. If
#'        more than one column were used to define the condition, this variable
#'        holds the effecive condition which is the combination of all applied
#'        single conditions. }
#' }
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
reduceReplicates <- function(file, select = "condition", datatrans = NULL) {

  # File format definiton
  fmtnames <- c("name", "time",  "value", "condition")
  fmtnamesnumber <- length(fmtnames)

  # Read data and sanity checks
  data <- read.csv(file)
  if (length(intersect(names(data), fmtnames)) != fmtnamesnumber) {
    stop(paste("Mandatory column names are:", paste(fmtnames, collapse = ", ")))
  }

  # Transform data if requested
  if (is.character(datatrans)) {
    x <- data$value
    data$value <- eval(parse(text = datatrans))
  }

  # Experiments are usually repeated multiple times possibly under different
  # conditions. The column "Condition" in the data thus groups the data per
  # condition. However, sometimes, a more fine grained grouping is desirable. In
  # this case, any number of additional columns can be append to the data. These
  # columns are referred to as "condition identifier". Which of the condition
  # identifiers are used to do the grouping is user defined by giving their
  # names in <select>. The mandatory column "Condition" is always used. The
  # total set of different conditions is thus defined by all combinations of
  # values occuring in the condition identifiers named for grouping. Mean and
  # variance is computed for each condition by averaging over measurements
  # recorded at the same time point. New conditions names are derived by merging
  # all conditions which were used in mean and std.
  select <- unique(c("name", "time", "condition", select))
  condidnt <- Reduce(paste, subset(data, select = select))
  conditions <- unique(condidnt)
  reduct <- do.call(rbind, lapply(conditions, function(cond) {
    conddata <- data[condidnt == cond,]
    mergecond <- Reduce(paste, conddata[1, setdiff(select, c("name", "time"))])
    data.frame(time = conddata[1, "time"],
               value = mean(conddata[, "value"]),
               sigma = sd(conddata[, "value"])/sqrt(nrow(conddata)),
               n = nrow(conddata),
               name = conddata[1, "name"],
               condition = mergecond)
  }))

  return(reduct)
}



#' Fit an error model
#'
#' @description Fit an error model to reduced replicate data, see
#'   \code{\link{reduceReplicates}}.
#'
#' @param data Reduced replicate data, see \code{\link{reduceReplicates}}. Need 
#'   columns "value", "sigma", "n".
#' @param factors \option{data} is pooled with respect to the columns named
#'   here, see Details.
#' @param errorModel Character vector defining the error model in terms of the variance. 
#'   Use \kbd{x} to reference the independend variable, see Details.
#' @param par Inital values for the parameters of the error model.
#' @param plotting If TRUE, a plot of the pooled variance together with the fit
#'   of the error model is shown.
#' @param blather If TRUE, additional information is returned, such as fit parameters 
#'  and sigmaLS (original sigma given in input data).
#' @param ... Parameters handed to the optimizer \code{\link{optim}}.
#'
#' @details The variance estimator using \eqn{n-1} data points is \eqn{chi^2}
#'   distributed with \eqn{n-1} degrees of freedom. Given replicates for
#'   consecutive time points, the sample variance can be assumed a function of
#'   the sample mean. By defining an error model which must hold for all time
#'   points, a maximum likelihood estimator for the parameters of the error
#'   model can be derived. The parameter \option{errorModel} takes the error
#'   model as a character vector, where the mean (independent variable) is
#'   refered to as \kbd{x}.
#'
#'   It is desireable to estimate the variance from many replicates. The
#'   parameter \option{data} must provide one or more columns which define the
#'   pooling of data. In case more than one column is announced by
#'   \option{factors}, all combinations are constructed. If, e.g.,
#'   \option{factors = c("condition", "name")} is used, where "condition" is
#'   "a", "b", "c" and repeating and "name" is "d", "e" and repeating, the
#'   effective conditions used for pooling are "a d", "b e", "c d", "a e", "b
#'   d", and "c e".
#'
#'   By default, a plot of the pooled data, sigma and its confidence bound at
#'   68\% and 95\% is shown.
#'
#' @return Returned by default is a data frame with columns as in \option{data}, 
#'   but with the sigma values replaced by the derived values, obtained by evaluating 
#'   the error model with the fit parameters. 
#'   
#'   If the blather = TRUE option is chosen, fit values of the parameters of the error
#'   model are appended, with the column names equal to the parameter names. 
#'   The error model is appended as the attribute "errorModel".
#'   Confidence bounds for sigma at confidence level 68\% and 95\% are
#'   calculated, their values come next in the returned data frame. Finally, the
#'   effective conditions are appended to easily check how the pooling was done.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
#' @importFrom stats D approx optim qchisq sd time
fitErrorModel <- function(data, factors, errorModel = "exp(s0)+exp(srel)*x^2",
                          par = c(s0 = 1, srel = .1), plotting = TRUE, blather = FALSE, ...) {

  # Assemble conditions
  condidnt <- Reduce(paste, subset(data, select = factors))
  conditions <- unique(condidnt)


  # Fit error model
  nColData <- ncol(data)
  dataErrorModel <- cbind(data, as.list(par))

  for (cond in conditions) {
    subdata <- dataErrorModel[condidnt == cond,]
    x <- subdata$value
    n <- subdata$n
    y <- subdata$sigma*sqrt(n)

    obj <- function(par) {
      value <- with(as.list(par), {
        z <- eval(parse(text = errorModel))
        sum(log(z)-log(dchisq((n-1)*(y^2)/z, df = n-1)), na.rm = TRUE)
      })
      return(value)
    }

    fit <- optim(par = par, fn = obj, ...)
    sigma <- sqrt(with(as.list(fit$par), eval(parse(text = errorModel))))
    dataErrorModel[condidnt == cond, ]$sigma <- sigma 
    dataErrorModel[condidnt == cond, -(nColData:1)] <- data.frame(as.list(fit$par))
  }


  # Calculate confidence bounds about sigma
  p68 <- (1-.683)/2
  p95 <- (1-.955)/2
  dataErrorModel$cbLower68 <- dataErrorModel$sigma^2*qchisq(p = p68, df = dataErrorModel$n-1)/(dataErrorModel$n-1)
  dataErrorModel$cbUpper68 <- dataErrorModel$sigma^2*qchisq(p = p68, df = dataErrorModel$n-1, lower.tail = FALSE)/(dataErrorModel$n-1)
  dataErrorModel$cbLower95 <- dataErrorModel$sigma^2*qchisq(p = p95, df = dataErrorModel$n-1)/(dataErrorModel$n-1)
  dataErrorModel$cbUpper95 <- dataErrorModel$sigma^2*qchisq(p = p95, df = dataErrorModel$n-1, lower.tail = FALSE)/(dataErrorModel$n-1)


  # Assemble result
  dataErrorModel <- cbind(dataErrorModel, condidnt, sigmaLS = data$sigma)
  attr(dataErrorModel, "errorModel") <- errorModel


  # Plot if requested
  if (plotting) {
    print(ggplot(dataErrorModel, aes(x=value)) +
            geom_point(aes(y=sigmaLS^2*(n))) +
            geom_line(aes(y=sigma^2)) +
            geom_ribbon(aes(ymin=cbLower95, ymax=cbUpper95), alpha=.3) +
            geom_ribbon(aes(ymin=cbLower68, ymax=cbUpper68), alpha=.3) +
            ylab("variance") +
            facet_wrap(~condidnt, scales = "free") +
            scale_y_log10() +
            theme_dMod()
    )}
  
  # Return standard error of the mean
  dataErrorModel$sigma <- dataErrorModel$sigma/sqrt(dataErrorModel$n)
  data$sigma <- dataErrorModel$sigma
  if(blather)
    return(dataErrorModel)
  else 
    return(data)
}



#' Concatenate parameter lists
#'
#' @description Fitlists carry an fit index which must be held unique on merging
#' multiple fitlists.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @rdname parlist
#' @export
#' @export c.parlist
c.parlist <- function(...) {
    m_fits <- lapply(list(...), unclass)
    m_fits <- do.call(c, m_fits)
    m_parlist <- mapply(function(fit, idx) {
      if (is.list(fit)) fit$index <- idx
      return(fit)
      }, fit = m_fits, idx = seq_along(m_fits))
    
    return(as.parlist(m_parlist))
  }



#' Check which Python versions are installed on the system
#' 
#' @param version NULL or character. Check for specific version
#' @return Character vector with the python versions and where they are located.
#' @export
python.version.sys <- function(version = NULL) {
  
  # Which python versions are installed on the system
  m_sysPath <- strsplit(Sys.getenv("PATH"), ":")
  m_sysPath <- m_sysPath[[1]]
  m_python <- do.call(rbind, lapply(m_sysPath, function(p) {
    m_py <- dir(p, pattern = "^python[0-9,.]+$")
    if (length(m_py) != 0) {
      return(file.path(p, m_py))
    } else {
      return(NULL)
    }
  }))
  
  m_version <- strsplit(m_python, "python")
  m_version <- lapply(m_version, function(p) {
    return(p[2])
  })
  
  m_python <- as.data.frame(m_python)
  names(m_python) <- m_version
  
  
  if (is.null(version)) {
    return(m_python)  
  } else {
    # Is requested version available
    if (any(m_version == version)) {
      return(as.character(m_python[[version]]))
      attr(out, "version") <- m_version
    } else {
      return(NULL)
    }
  }
}


#' Get the Python version to which rPython is linked
#' 
#' @return The Python version and additional information
#' @export
python.version.rpython <- function() {
  rPython::python.exec(c("def ver():", "\timport sys; return list(sys.version_info)"))
  m_info <- as.data.frame(rPython::python.call("ver"))
  names(m_info) <- c("major", "minor", "micro", "releselevel", "serial")
  
  m_version <- paste0(m_info[[1]], ".", m_info[[2]])
  attr(m_version, "info") <- m_info
  
  return(m_version)
}


#' Check if rPython comes with the correct Python version
#' 
#' @description rPython is liked against a certain Python version found on the system.
#' If Python code called from R requires a specific Python version, the rPython package
#' needs to be reinstalled. This functions helps to do this in one line.
#' 
#' @param version character indicating the requested Python version
#' 
#' @return TRUE if rPython is linked against the requested version. Otherwise, the user
#' is asked if rPython should be reinstalled with the correctly linked Python version.
#' @export
python.version.request <- function(version) {
  
  # Is rPythen installed and linked against requested python version?
  m_installed <- "rPython" %in% installed.packages()[, 1]
  if (m_installed) {
    m_curVersion <- python.version.rpython()
    if (m_curVersion == version)
      return(TRUE)
  }
  
  # rPython not installed or linked agains wrong python version.
  m_sysVersion <- python.version.sys(version)
  if (is.null(m_sysVersion)) {
    msg <- paste0("Requested python version ", version, " not available\n",
                  "Your options are\n")
    cat(msg)
    print(python.version.sys())
  } else {
    # If rPython is already installed, double check with user
    if (m_installed) {
      msg <- paste0("rPython is installed on your system using python ", m_curVersion, "\n",
                    "Proceeding the installation enables ", version, " for R\n",
                    "This will prevent python programms needing version ", m_curVersion, " from executing\n",
                    "Do you want to abort [1] or procede [2] with the installation?\n")
      cat(msg)
      m_go <- scan(nmax = 1, what = integer(), quiet = TRUE)
      if (m_go != 2) {
        stop("Installation aborted")
      }
    }
    
    try(detach("package:rPython", unload = TRUE), silent = TRUE)
    Sys.setenv(RPYTHON_PYTHON_VERSION = version)
    install.packages("rPython")
  }
}