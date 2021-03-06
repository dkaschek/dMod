% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{confint.parframe}
\alias{confint.parframe}
\title{Profile uncertainty extraction}
\usage{
\method{confint}{parframe}(object, parm = NULL, level = 0.95, ..., val.column = "data")
}
\arguments{
\item{object}{object of class \code{parframe}, returned from \link{profile} function.}

\item{parm}{a specification of which parameters are to be given confidence intervals, 
either a vector of numbers or a vector of names. If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{...}{not used right now.}

\item{val.column}{the value column used in the parframe, usually 'data'.}
}
\description{
extract parameter uncertainties from profiles
}
