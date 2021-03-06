% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steadyStates.R
\name{python_version_request}
\alias{python_version_request}
\title{Check if rPython comes with the correct Python version}
\usage{
python_version_request(version)
}
\arguments{
\item{version}{character indicating the requested Python version}
}
\value{
TRUE if rPython is linked against the requested version. Otherwise, the user
is asked if rPython should be reinstalled with the correctly linked Python version.
}
\description{
rPython is liked against a certain Python version found on the system.
If Python code called from R requires a specific Python version, the rPython package
needs to be reinstalled. This functions helps to do this in one line.
}
