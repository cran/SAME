#'
#' A Time-to-event dataset containing the time and other attributes of 643 patients.
#'
#' @format A data frame with 643 rows and 6 variables:
#' \describe{
#'   \item{T1}{binary variable, receive treatment 1=1, not receive treatment 1=0}
#'   \item{T2}{binary variable, receive treatment 2=1, not receive treatment 2=0}
#'   \item{B1}{binary variable, biomarker 1 positive=1, biomarker 1 negative=0}
#'   \item{B2}{binary variable, biomarker 2 positive=1, biomarker 2 negative=0}
#'   \item{death}{the status indicator, alive=0, dead=1}
#'   \item{surv}{survival time or follow up time}
#'   \item{survtime}{survival time or follow up time}
#'   \item{treatments}{categorical vairable, indicating treatments received}
#'   ...
#' }
"example.2"
