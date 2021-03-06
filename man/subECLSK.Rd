% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subECLSK.R
\docType{data}
\name{subECLSK}
\alias{subECLSK}
\title{A subset of the ECLSK (Early Childhood Longitudinal Study – Kindergarten) cohort data.}
\format{A data frame with 600 rows (participants) and 11 columns (variables).}
\usage{
data(subECLSK)
}
\description{
The dataset contains the reading and mathematics scores and other attributes of
600 kindergarteners. The variables are as follows.
}
\details{
\itemize{
\item dobmm. The date of birth measured in months.
\item gender. 1=Male; 2=Female
\item race.
\itemize{
\item{1 = White,non-Hispanic}
\item{2 = Black or African American, non-Hispanic}
\item{3 = Hispanic, race specified}
\item{4 = Hispanic, race no specified}
\item{5 = Asian}
\item{6 = Native Hawaiin, other Pacific islander}
\item{7 = American Indian or Alaska native}
\item{8 = More than one race, non-Hispanic}
}
\item readingIRT. Reading IRT (Item Response Theory) scaled score.
\item mathIRT. Mathematics IRT (Item Response Theory) scaled score.
\item numsib. Number of siblings in household.
\item parentedu. Parent highest education level.
\item ses. Continuous socioeconomic status measure.
\item relAge. The relative age of the participant entering kindergarten measured in months.
\item PredEnt. The predicted relative kindergarten entrance age.
\item mis.ind.read. Is the reading IRT score missing? 1=Yes, 0=No.
\item mis.ind.math. Is the math IRT score missing? 1=Yes, 0=No.
}
}
\examples{
\donttest{
data(subECLSK)
}

}
\references{
Zhong, S. Y., & Hoxby, C. M. (2012).
The effects of relative age on early childhood academic achievement:
how they differ between gender and change across time.
(Unpublished Honors Thesis) Stanford University, Stanford, CA.
\href{http://economics. stanford. edu/_les/StephanieYuechenZhongHonorsThesis2012.}{http://economics. stanford. edu/_les/StephanieYuechenZhongHonorsThesis2012.}

Tourangeau, K., Nord, C., Lê, T., Sorongon, A. G., & Najarian, M. (2009).
Early childhood longitudinal study, kindergarten class of 1998–99 (ECLS-K),
combined User's manual for the ECLS-K eighth-grade and K–8 full sample data files
and electronic codebooks (NCES2009–004).
National Center for Education Statistics, Institute of Education Sciences,
U.S. Departmentof Education. Washington DC.
}
\keyword{datasets}
