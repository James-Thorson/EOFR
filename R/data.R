
#' Data to demonstrate EOFR package
#'
#' Data sufficient to demonstrate a model jointly applying empirical orthogonal function to physical variables
#'   and a linear regression to a biological (or other) time-series, where modes of variability from physical
#'   variables are used as predictors in the time-series model
#'
#' \itemize{
#'   \item physical_data data-frame of biological sampling data
#'   \item stock_recruit_data stock-recruit records for Pacific cod, with columns for log-recruits-per-spawner (response), spawning biomass (predictor), and year (associated with year in physical variable), where recruits-per-spawning biomass is lagged to associate spawning biomass with estimated age-0 abundance in each year
#' }
#'
#' @name EOFR_example
#' @docType data
#' @usage data(EOFR_example)
#' @keywords data
NULL

