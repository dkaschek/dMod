% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normIndiv.R
\name{getParGrids}
\alias{getParGrids}
\title{Get est.grid and fixed.grid}
\usage{
getParGrids(mytrafo, mytrafoL, mycondition.grid, SS_pars = NULL)
}
\arguments{
\item{mytrafo}{base trafo}

\item{mytrafoL}{condition specific branched trafo list}

\item{mycondition.grid}{condition.grid with condition names as rownames e.g. as output from attr(datalist, "condition.grid")}

\item{SS_pars}{parameters determined by the steady state}
}
\description{
Get est.grid and fixed.grid
}
\author{
Svenja Kemmer
}
