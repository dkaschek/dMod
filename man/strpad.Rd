% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strings.R
\name{strpad}
\alias{strpad}
\title{Pad string to desired width}
\usage{
strpad(string, width, where = "right", padding = " ", autoelide = FALSE)
}
\arguments{
\item{string}{String to pad}

\item{width}{Desired width of padded string}

\item{where}{Padding can be inserted to the right or left of <string>.
Default to 'right'.}

\item{padding}{A single character with with the padding space is filled.
Defaults to blank ' ' yielding invisible padding.}

\item{autoelide}{If TRUE, <string> is elided if it is wider than <width>. The
position of eliding follows <where>. Defaults to FALSE.}
}
\value{
Padded string of length <width>.
}
\description{
Pad string to desired width
}
\author{
Wolfgang Mader, \email{Wolfgang.Mader@fdm.uni-freiburg.de}
}
