% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotFluxes}
\alias{plotFluxes}
\title{Plot Fluxes given a list of flux Equations}
\usage{
plotFluxes(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...)
}
\arguments{
\item{pouter}{parameters}

\item{x}{The model prediction function \code{x(times, pouter, fixed, ...)}}

\item{times}{Numeric vector of time points for the model prediction}

\item{fluxEquations}{list of chars containing expressions for the fluxes, 
if names are given, they are shown in the legend. Easy to obtain via \link{subset.eqnlist}, see Examples.}

\item{nameFlux}{character, name of the legend.}

\item{...}{Further arguments going to x, such as \code{fixed} or \code{conditions}}
}
\value{
A plot object of class \code{ggplot}.
}
\description{
Plot Fluxes given a list of flux Equations
}
\examples{
\dontrun{

plotFluxes(bestfit, x, times, subset(f, "B"\%in\%Product)$rates, nameFlux = "B production")
}
}
