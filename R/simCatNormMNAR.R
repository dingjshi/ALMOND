#' A simulated dataset 2.
#'
#' A simulated dataset with the categorical treatment variable. The outcome variable
#' is normally-distributed and is missing not at random (MNAR, e.g., dropout or attrition).
#' The dataset contains a data frame with 600 rows (participants) and 4 columns (variables).
#' The variables are as follows.
#'
#' \itemize{
#' \item outcome. The hypothetical causal outcome variable.
#' \item treatment. The hypothetical causal treatment variable.
#' \item instrument. The hypothetical instrumental variable.
#' \item mis.ind. Is the outcome variable value missing? 1=Yes, 0=No.
#' }
#'
#' @name simCatNormMNAR
#'
#' @docType data
#'
#' @usage data(simCatNormMNAR)
#'
#' @format A data frame with 600 rows and 4 columns .
#'
#' @keywords datasets
#'
#' @examples
#' \donttest{
#' data(simCatNormMNAR)
#' }
#'
NULL
#------
