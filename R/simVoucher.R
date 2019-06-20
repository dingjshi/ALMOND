#' Simulated data of a public voucher program.
#'
#' The dataset is simulated from a study of the effect of a public voucher program
#' (Currie and Yelowitz, 2000). The dataset contains the status of the public voucher program
#' participation and other attributes. The variables are as follows.
#'
#' \itemize{
#' \item nbhoodRating. Rating of neighborhood.
#' \item hmRating. Rating of home.
#' \item voucher. Participated in the voucher program? 1=Yes; 0=No.

#' \item extraBed. The household has an extra bedroom? 1=Yes; 0=No.
#' \item girl. Is the child a girl? 1=Yes; 0=No.
#' \item childAge. Age of the child in the household.
#' \item headFemale. Is the household head a female? 1=Yes; 0=No.
#' \item headAge. Age of the household head.
#' \item headMarried. Is the household head married? 1=Yes; 0=No.
#' \item hdEd. Years of education of the household head.
#' \itemize{
#' \item{1 = Between 9 and 11.}
#' \item{2 = Equals 12.}
#' \item{3 = Between 13 and 15.}
#' \item{4 = Greater than 16.}
#' }
#'
#' @name simVoucher
#'
#' @docType data
#'
#' @usage data(simVoucher)
#'
#' @format A data frame with 1954 rows (participants) and 10 columns (variables).
#'
#' @keywords datasets
#'
#' @references Currie, J., & Yelowitz, A. (2000).
#' Are public housing projects good for kids? \emph{Journal of public economics},
#' \emph{75}(1), 99-124.
#'
#' @examples
#' \donttest{
#' data(simVoucher)
#' }
#'
NULL
#------
