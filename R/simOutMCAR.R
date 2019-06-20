#' A simulated dataset 7.
#'
#' A simulated dataset with the continuous treatment variables. The outcome variable contains
#' outliers and is missing completely at random (MCAR).
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
#' @name simOutMCAR
#'
#' @docType data
#'
#' @usage data(simOutMCAR)
#'
#' @format A data frame with 600 rows and 4 columns .
#'
#' @keywords datasets
#'
#' @examples
#' \donttest{
#' data(simOutMCAR)
#' }
#'
NULL
#------
