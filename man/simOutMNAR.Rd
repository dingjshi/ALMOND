% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simOutMNAR.R
\docType{data}
\name{simOutMNAR}
\alias{simOutMNAR}
\title{A simulated dataset 8.}
\format{A data frame with 600 rows and 4 columns .}
\usage{
data(simOutMNAR)
}
\description{
A simulated dataset with the continuous treatment variable. The outcome variable contains
outliers and is missing not at random (MNAR, e.g., dropout or attrition).
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
data(simOutMNAR)
}

}
\keyword{datasets}
