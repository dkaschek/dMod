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
#'
#' @export
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
#'
#' @export
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



#' Select parameter values with lowest Chi^2 among profiles.
#'
#' @param prf A profiles as returned from \link[dMod]{profile}.
#' @param context If TRUE, the chi^2 and other context of the profile is
#'                returned. This parameter defaults to FALSE in which case the
#'                output can be used as outer parameters directly.
#'
#' @return Parameter values with lowest chi^2 w/ or w/o profile context. If all
#'         profiles are invalid, NULL is returned.
#'
#' @details
#' On profiling a model, parameter values yielding a lower chi^2 than the one of
#' the fit providing the optimal values for the profiles might be encountered.
#' This function extracts the set of parameter values possessing the lowest
#' chi^2 among all profiles. Profiles which are invalid due to e.g., integration
#' problem on calculating them, are ignored. If all profiles are invalid, NULL
#' is returned.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export

plSelectMin <- function(prf, context = FALSE) {

  # Remove invalid profiles.
  apt <- sapply(prf, class) == "matrix"
  if (!any(apt)) {
    return(NULL)
  } else {
    prf <- prf[apt]
  }

  # Minium chi^2 parameter values sets per profile.
  chi2MinAllProfiles <- sapply(prf, function(species) {
    return(species[which.min(species[, "value"]), ])
  })

  # Minimum chi^2 parameter values set accross all profiles.
  chi2MinBest <- chi2MinAllProfiles[, which.min(chi2MinAllProfiles[1, ])]

  if (context) {
    return(chi2MinBest)
  } else {
    return(chi2MinBest[-(1:4)])
  }

}



#' Non-Linear Optimization, multi start
#'
#' @description
#' Wrapper around \code{\link{trust}} allowing for multiple fits from randomly
#' chosen initial values.
#'
#' @param objfun Objective function, see \code{\link{trust}}.
#' @param center Parameter values around which the initial values for each fit
#'        are randomly sampled. The initial values handed to \link{trust} are
#'        the sum of center and the output of <samplefun>, center + <samplefun>.
#'        See \code{\link{trust}}, parinit.
#'
#' @param rinit Starting trust region radius, see \code{\link{trust}}.
#' @param rmax Maximum allowed trust region radius, see \code{\link{trust}}.
#' @param fits Number of fits (jobs).
#' @param cores Number of cores for job parallelization.
#' @param samplefun Function to sample random initial values. It is assumed,
#'        that <samplefun> has a named parameter "n" which defines how many
#'        random numbers are to be returned, such as for \code{\link{rnorm}},
#'        which is also the default sampling function.
#' @param logfile Name of the file to which all jobs log their output. The file
#'        is handy to investigate the different jobs in some detail. Since the
#'        jobs are carried out in parallel, their output may occurre in
#'        non-consecutive order. At the end of the file, a summary of the fits
#'        is given.
#' @param fitsfile Name of the file to which the result of all completed fits
#'        are written to. An empy string "" suppresses the write.
#' @param stats If true, the same summary statistic as written to the logfile is
#'        printed to command line.
#' @param msgtag A string prepending the logging output written to file.
#' @param ... Additional parameters which are handed to trust() or samplefun()
#'        by matching parameter names. All remaining parameters are handed to
#'        the objective function objfun().
#'
#' @details
#' By running multiple fits starting with randomly chosen inital parameters, the
#' chisquare landscape can be explored using a deterministic optimizer. In this
#' case, \code{\link{trust}} is used for optimization. The standard procedure to
#' obtain random initial values is to sample random variables from a uniform
#' distribution (\code{\link{rnorm}}) and adding these to <center>. It is,
#' however, possible, to employ any other sampling strategy by handing the
#' respective function to mstrust(), <samplefun>.
#'
#' All started fits are either faulty, aborted, or complete. Faulty fits return
#' a "try-error" object and fail somewhere outside trust(). Aborted fits fail
#' withing trust(), and complete fits return from trust() correctly. Completed
#' fits can still be unconverged, in case the maximum number of iteration is
#' reached before the convergence criterion.
#'
#' @return fitlist A data frame of all completed fits sorted by their chi^2.
#'
#' @seealso
#' The optimization is carried out by \code{\link{trust}}.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export
mstrust <- function(objfun, center, rinit = .1, rmax = 10, fits = 20, cores = 1,
                    samplefun = "rnorm", logfile = "mstrust.log",
                    fitsfile = "fitlist.rda", stats = FALSE, msgtag = "",
                    ...) {

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
  # Second, check what trust() and samplefun() accept and check for names clash.
  # Third, whatever is unused is passed to the objective function objfun().
  nameslocal <- c("center", "fits", "cores", "samplefun", "logfile", "msgtag", "stats", "writeres")
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
  for (name in namestrust){
    argstrust[[name]] <- argslist[[name]]
  }


  # Apply trust optimizer in parallel
  # The error checking leverages that mclappy runs each job in a try().
  file.create(argslist$logfile) #Truncate log file
  logfile <- file(argslist$logfile, open = "a")

  fitlist <- mclapply(1:fits, function(i) {
    argstrust$parinit <- center + do.call(samplefun, argssample)
    fit <- do.call(trust, argstrust)
    fit$index = i

    # Reporting
    # With concurent jobs and everyone reporting, this is a classic race
    # condition. Assembling the message beforhand lowers the risk of interleaved
    # output to the log.
    msgTag <- argslist$msgtag
    msgSep <- "-------"
    if (any(names(fit) == "error")) {
      msg <- paste0(msgTag, msgSep, "\n",
                    msgTag, "Fit ", i, " failed after ", fit$iterations, " iterations with error\n",
                    msgTag, "--> ", fit$error,
                    msgTag, msgSep, "\n")

      writeLines(msg, logfile)
      flush(logfile)
    } else {
      msg <- paste0(msgTag, msgSep, "\n",
                    msgTag, "Fit ", i, " completed\n",
                    msgTag, "--> iterations : ", fit$iterations, "\n",
                    msgTag, "-->  converged : ", fit$converged, "\n",
                    msgTag, "-->      chi^2 : ", round(fit$value, digits = 2), "\n",
                    msgTag, msgSep)

      writeLines(msg, logfile)
      flush(logfile)
    }

    return(fit)
  },mc.preschedule=FALSE, mc.silent = FALSE, mc.cores=cores)
  close(logfile)


  # Cull failed, aborted, and completed fits
  # Failed jobs return an object of class "try-error". The reason for these
  # failures are unknown to me.
  # Aborted fits, in contrast, return a list of results from trust(), where one
  # name of the list is error holding an object of class "try-error". These
  # abortions are due to errors which are captured within trust().
  # Completed fits return with a valid result list from trust(), with error is
  # not part of its names. These fits, however, can still be unconverged, if the
  # maximim number of iterations was the reason for the return of trust().
  # Failed: try-error, find and remove
  idxerr <- sapply(fitlist, function(fit) inherits(fit, "try-error"))
  fitlist <- fitlist[!idxerr]

  # Aborted: Due to some error condition handled inside trust()
  idxabrt <- sapply(fitlist, function(fit) {
    if (any(names(fit) == "error")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  # Properly returned from trust(), possibly unconverged.
  idxcmp <- !idxabrt


  # Stash results of completed fits in a data.frame
  complist <- lapply(fitlist[idxcmp], function(fit) {
    data.frame(
      index = fit$index,
      chisquare = fit$value,
      converged = fit$converged,
      iterations = fit$iterations,
      as.data.frame(as.list(fit$argument))
    )
  })
  compframe <- do.call(rbind, complist)
  if (nrow(compframe) > 0) {
    compframe <- compframe[order(compframe$chisquare),]
  }


  # Wrap up
  # Write out results
  if (nchar(fitsfile) > 0) {
    save(fitlist, file = fitsfile)
  }

  # Show summary
  msg <- paste0("Mutli start trust summary\n",
                "Outcome   : Occurrence\n",
                "Faulty    :", sum(idxerr), "\n",
                "Aborted   :", sum(idxabrt), "\n",
                "Completed :", sum(idxcmp), "\n",
                "           -----------\n",
                "Total     :", sum(idxerr) + sum(idxabrt) + sum(idxcmp), paste0("[", fits, "]"), "\n")
  logfile <- file(argslist$logfile, open = "a")
  writeLines(msg, logfile)
  flush(logfile)
  close(logfile)

  if (stats) {
    cat(msg)
#     cat("Mutli start trust summary\n")
#     cat("Outcome    : Occurrence\n")
#     cat(" Faulty    :", sum(idxerr), "\n")
#     cat(" Aborted   :", sum(idxabrt), "\n")
#     cat(" Completed :", sum(idxcmp), "\n")
#     cat("             -----------\n")
#     cat(" Total     :", sum(idxerr) + sum(idxabrt) + sum(idxcmp), paste0("[", fits, "]"), "\n")
  }

  return(compframe)
}



#' Select best fit.
#'
#' @description
#' Select the fit with lowest chi^2 form the result of \code{\link{mstrust}}.
#'
#' @param fitlist A data frame of fits as returned from \code{\link{mstrust}}.
#'        The data frame does not need to be ordered and can include unconverged
#'        fits.
#'
#' @return The converged fit with lowest chisquare as a named numeric vector.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @export

msbest <- function(fitlist) {
  fitlistconv <- fitlist[fitlist$converged == TRUE,]
  if (is.null(nrow(fitlistconv))) {
    return(NULL)
  }

  idxbest <- order(fitlistconv$chisquare)
  best <- fitlistconv[idxbest[1],]
  best <- unlist(best[1, -(4:1)])

  return(best)
}
