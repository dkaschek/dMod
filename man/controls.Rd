% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{controls}
\alias{controls}
\alias{controls.objfn}
\alias{controls.fn}
\alias{controls<-}
\alias{controls<-.objfn}
\alias{controls<-.fn}
\title{List, get and set controls for different functions}
\usage{
controls(x, ...)

\method{controls}{objfn}(x, name = NULL, ...)

\method{controls}{fn}(x, condition = NULL, name = NULL, ...)

controls(x, ...) <- value

\method{controls}{objfn}(x, name, ...) <- value

\method{controls}{fn}(x, condition = NULL, name, ...) <- value
}
\arguments{
\item{x}{function}

\item{...}{arguments going to the appropriate S3 methods}

\item{name}{character, the name of the control}

\item{condition}{character, the condition name}

\item{value}{the new value}
}
\value{
Either a print-out or the values of the control.
}
\description{
Applies to objects of class \code{objfn},
\code{parfn}, \code{prdfn} and \code{obsfn}. Allows to manipulate
different arguments that have been set when creating the
objects.
}
\details{
If called without further arguments, \code{controls(x)} lists the
available controls within an object. Calling \code{controls()} with \code{name}
and \code{condition} returns the control value. The value can be overwritten. If
a list or data.frame ist returned, elements of those can be manipulated by the
\code{$}- or \code{[]}-operator.
}
\examples{
## parfn with condition
p <- P(eqnvec(x = "-a*x"), method = "implicit", condition = "C1")
controls(p)
controls(p, "C1", "keep.root")
controls(p, "C1", "keep.root") <- FALSE

## obsfn with NULL condition
g <- Y(g = eqnvec(y = "s*x"), f = NULL, states = "x", parameters = "s")
controls(g)
controls(g, NULL, "attach.input")
controls(g, NULL, "attach.input") <- FALSE
}
