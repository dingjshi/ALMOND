% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simNormMCAR.R
\docType{data}
\name{simNormMCAR}
\alias{simNormMCAR}
\title{A simulated dataset 5.}
\format{A data frame with 600 rows and 4 columns .}
\usage{
data(simNormMCAR)
}
\description{
A simulated dataset with the continous treatment variable. The outcome variable is
normally-distributed and is missing completely at random (MCAR).
The dataset contains a data frame with 600 rows (participants) and 4 columns (variables).
The variables are as follows.
}
\details{
\itemize{
\item outcome. The hypothetical causal outcome variable.
\item treatment. The hypothetical causal treatment variable.
\item instrument. The hypothetical instrumental variable.
\item mis.ind. Is the outcome variable value missing? 1=Yes, 0=No.
}
}
\examples{
\donttest{
data(simNormMCAR)
}

}
\keyword{datasets}
