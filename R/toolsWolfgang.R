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
