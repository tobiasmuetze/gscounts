% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{data}
\name{hospitalizations}
\alias{hospitalizations}
\title{Hospitalizations}
\format{
A data frame with 2323 rows and 4 variables
}
\usage{
data(hospitalizations)
}
\description{
A dataset containing the hospitalization times of 1980 patients:
}
\details{
\itemize{
  \item treatment. Treatment identifier. 
  \item pat. Patient identifier. Unique within \code{treatment}
  \item t_recruit. Recruitment time of patient into the clinical trial.
  \item eventtime. Event time of hospitalization. \code{NA} corresponds to no event.
}
}
\keyword{datasets}
