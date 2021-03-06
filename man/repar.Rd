% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{repar}
\alias{repar}
\title{Reparameterization}
\usage{
repar(expr, trafo = NULL, ..., reset = FALSE)
}
\arguments{
\item{expr}{character of the form \code{"lhs ~ rhs"} where \code{rhs}
reparameterizes \code{lhs}. Both \code{lhs} and \code{rhs}
can contain a number of symbols whose values need to be passed by the \code{...} argument.}

\item{trafo}{character or equation vector or list thereof. The object where the replacement takes place in}

\item{...}{pass symbols as named arguments}

\item{reset}{logical. If true, the trafo element corresponding to lhs is reset according to rhs. 
If false, lhs wherever it occurs in the rhs of trafo is replaced by rhs of the formula.}
}
\value{
an equation vector with the reparameterization.
}
\description{
Reparameterization
}
\details{
Left and right-hand side of \code{expr} are searched for symbols. If separated by
"_", symbols are recognized as such, e.g. in \code{Delta_x} where the symbols are 
"Delta" and "x". Each symbol for which values (character or numbers) are passed by the
\code{...} argument is replaced.
}
\examples{
innerpars <- letters[1:3]
constraints <- c(a = "b + c")
mycondition <- "cond1"

trafo <- repar("x ~ x", x = innerpars)
trafo <- repar("x ~ y", trafo, x = names(constraints), y = constraints)
trafo <- repar("x ~ exp(x)", trafo, x = innerpars)
trafo <- repar("x ~ x + Delta_x_condition", trafo, x = innerpars, condition = mycondition)
}
