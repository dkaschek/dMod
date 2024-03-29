% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolsSvenja.R
\name{plotArray}
\alias{plotArray}
\title{Plot an array of trajectories along the profile of a parameter}
\usage{
plotArray(
  par,
  profs,
  prd,
  times,
  direction = c("up", "down"),
  covtable,
  ...,
  nsimus = 4
)
}
\arguments{
\item{par}{Character of parameter name for which the array should be generated.}

\item{profs}{Lists of profiles as being returned by \link{profile}.}

\item{prd}{Named list of matrices or data.frames, usually the output of a prediction function
as generated by \link{Xs}.}

\item{times}{Numeric vector of time points for the model prediction.}

\item{direction}{Character "up" or "down" indicating the direction the value should be traced along the profile starting at the bestfit value.}

\item{covtable}{Optional covariate table or condition.grid necessary if subsetting is required.}

\item{...}{Further arguments for subsetting the plot.}

\item{nsimus}{Number of trajectories/ simulation to be calculated.}
}
\value{
A plot object of class \code{ggplot}.
}
\description{
Plot an array of trajectories along the profile of a parameter
}
\examples{
\dontrun{
 plotArray("myparameter", myprofiles, g*x*p, seq(0, 250, 1), 
    "up", condition.grid, name == "ProteinA" & condition == "c1") 
}
}
\author{
Svenja Kemmer, \email{svenja.kemmer@fdm.uni-freiburg.de}
}
