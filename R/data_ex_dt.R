#' Example Survival Data for Cluster-Randomized Trials
#'
#' This dataset (`ex_dt`) contains simulated survival data for use in
#' doubly robust survival analysis in cluster-randomized trials.
#'
#' @format A data frame with multiple columns:
#' \describe{
#'   \item{time}{Survival time (numeric)}
#'   \item{Delta}{Event indicator (1-3 = multi-state event, 0 = censored)}
#'   \item{Z1}{A continuous individual-level covariate}
#'   \item{Z2}{A binary individual-level covariates}
#'   \item{W1}{A continuous cluster-level covariate}
#'   \item{W2}{A binary cluster-level covariates}
#'   \item{trt}{Treatment assignment (0 = control, 1 = treated)}
#'   \item{M}{Cluster ID}
#'   \item{id}{individual id}
#' }
#'
#' @source Simulated data for package examples.
#'
#' @examples
#' data(ex_dt)  # Load dataset
#' head(ex_dt)
#'
#' @usage data(ex_dt)
#' @keywords datasets
#' @name ex_dt
#' @docType data
NULL
