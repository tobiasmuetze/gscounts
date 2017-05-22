#' gscounts
#' 
#' Design and monitoring of group sequential designs with negative binomial data.
#' 
#' @docType package
#' @author Tobias Muetze <tobias.muetze@outlook.com>
#' @importFrom Rcpp evalCpp
#' @useDynLib gscounts
#' @name gscounts
NULL  


#' Hospitalizations
#'
#' A dataset containing the hospitalization times of 1980 patients:
#'
#' \itemize{
#'   \item treatment. Treatment identifier. 
#'   \item pat. Patient identifier. Unique within \code{treatment}
#'   \item t_recruit. Recruitment time of patient into the clinical trial.
#'   \item eventtime. Event time of hospitalization. \code{NA} corresponds to no event.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hospitalizations
#' @usage data(hospitalizations)
#' @format A data frame with 2323 rows and 4 variables
NULL