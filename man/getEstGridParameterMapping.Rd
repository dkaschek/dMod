% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normIndiv.R
\name{getEstGridParameterMapping}
\alias{getEstGridParameterMapping}
\title{Get Parameter mappings outer.estgrid - inner.estgrid}
\usage{
getEstGridParameterMapping(x)
}
\arguments{
\item{x}{est.grid}
}
\value{
named character
}
\description{
Get Parameter mappings outer.estgrid - inner.estgrid
}
\examples{
est.grid <- data.frame(ID = 1:2,
                       condition = c("A", "B"),
                       k1 = c("k1_A", "k1_B"),
                       k2 = c("k2_A", "k2_B"),
                       k3 = c("k3", NA),
                       k4 = c("k4", "k4"),
                       stringsAsFactors = FALSE)
getEstGridParameterMapping(est.grid)
}
\author{
Daniel Lill (daniel.lill@physik.uni-freiburg.de)
}
