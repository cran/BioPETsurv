#' Example dataset for package 'BioPETsurv'
#'
#' A dataset containing values of two biomarkers and survival outcomes of 1533 individuals.
#'
#' @docType data
#'
#' @usage data(SurvMarkers)
#'
#' @format A data frame with 1533 rows and 4 variables:
#' \describe{
#'   \item{time}{observed times of event or censoring}
#'   \item{event}{indicator of event; 0 means censored and 1 means event}
#'   \item{x1}{A modestly prognostic biomarker (concordance index=0.64)}
#'   \item{x2}{A strongly prognostic biomarker (concordance index=0.82)}
#' }
"SurvMarkers"

