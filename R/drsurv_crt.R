#' Doubly Robust Survival Analysis for Cluster-Randomized Trials
#'
#' This function estimates survival probabilities in cluster-randomized trials (CRTs)
#' using doubly robust estimation. It supports both restricted mean survival time (RMST)
#' and survival probability estimation. The function allows specifying marginal or frailty
#' models for the estimation.
#'
#' @param s_formula A formula specifying the survival outcome model (e.g., `Surv(time, event) ~ covariates`).
#' @param c_formula A formula specifying the censoring model (e.g., `Surv(time, event) ~ covariates + cluster(id)`).
#' @param trt A string indicating the treatment variable in `data`. The variable must be binary (`0/1`).
#' @param probs A numeric vector of length 2, specifying the treatment probabilities (`P(trt=1)`, `P(trt=0)`) in the CRT.
#' @param estimand A string indicating the estimand of interest. Options:
#'   - `"rmst"`: Restricted Mean Survival Time (default).
#'   - `"surv"`: Survival probability at specified time points (`res_time`).
#' @param scale A string specifying the scale of effect estimation:
#'   - `"RD"`: Risk Difference (default).
#'   - `"RR"`: Risk Ratio.
#' @param method A string indicating the estimation method:
#'   - `"marginal"`: Marginal Cox model (default).
#'   - `"frailty"`: Frailty Cox model.
#' @param res_time A numeric vector specifying the time points for estimation.
#' @param data A data frame containing the survival data for analysis.
#' @param id A string specifying the cluster/group ID variable in `data`.
#'
#' @return A list with two elements:
#'   - `s_c`: A data frame containing **cluster-level survival estimates** for each time point in `res_time`.
#'   - `s_i`: A data frame containing **individual-level survival estimates** for each time point in `res_time`.
#'
#' @import dplyr
#' @importFrom survival
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ex_dt)  # Data is loaded into `ex_dt`
#'
#' # Compute cluster sizes and treatment probabilities
#' sum_N <- ex_dt %>%
#'   group_by(M) %>%
#'   summarise(n())
#'
#' ex_dt$Ni <- rep(sum_N$`n()`, sum_N$`n()`)
#'
#' prop <- ex_dt %>%
#'   group_by(.data[["M"]]) %>%
#'   slice(1) %>%  # Extract first observation per cluster
#'   ungroup() %>%
#'   summarize(
#'     prop1 = mean(.data[["A"]] == 1),
#'     prop0 = mean(.data[["A"]] == 0)
#'   ) %>% as.data.frame()
#'
#' probs <- c(prop$prop1, prop$prop0)
#'
#' # Define survival and censoring models
#' s_formula <- Surv(obs_T,Delta) ~ Z1 + Z2 + W1 + W2 + Ni + cluster(M)
#' c_formula <- Surv(obs_T,Delta) ~ Z1 + Z2 + W1 + W2 + cluster(M)
#'
#' # Run doubly robust survival estimation
#' results <- drsurv_crt(
#'   s_formula = s_formula,
#'   c_formula = c_formula,
#'   trt = "A",
#'   probs = probs,
#'   estimand = "rmst",
#'   scale = "RD",
#'   method = "marginal",
#'   res_time = c(0.1, 0.5, 1, 1.5),  # Example time points
#'   data = ex_dt,
#'   id = "id"
#' )
#'
#' # View results
#' print(results$s_c)
#' print(results$s_i)
#' }
#'
#' @export
drsurv_crt <- function(s_formula, c_formula, trt, probs,
                       estimand="rmst", scale="RD", method="marginal", res_time,data,id){

  # Ensure formulas are actually formula objects
  if (!inherits(s_formula, "formula") || !inherits(c_formula, "formula")) {
    stop("Both s_formula and c_formula must be valid R formula objects.")
  }

  # Extract terms from formulas
  s_terms <- terms(s_formula)
  c_terms <- terms(c_formula)


  s_lhs <- all.vars(s_formula[[2]])  # Extracts Surv(time, event)
  c_lhs <- all.vars(c_formula[[2]])

  # Ensure Surv() format
  if (length(s_lhs) != 2 || length(c_lhs) != 2) {
    stop("Both formulas must be in the form Surv(time, event) ~ covariates.")
  }


  rel_s <- all.vars(s_formula)
  rel_c <- all.vars(c_formula)


  # Check if all variables exist in the dataset
  missing_vars <- setdiff(unique(c(rel_s,rel_c)), names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following variables are missing in the dataset:", paste(missing_vars, collapse = ", ")))
  }



  # Check if time is numeric
  if (!is.numeric(data[[s_lhs[1]]]) || !is.numeric(data[[c_lhs[1]]]) ) {
    stop(paste("The time variable must be numeric"))
  }

  if(s_lhs[1] != c_lhs[1]){
    stop(paste("The time variables for outocme and censoring should be the same"))

  }


  # Check if cluster() exists in c_formula
  cluster_var_c <- NULL
  cluster_var_s <- NULL
  if (any(grepl("cluster\\(", labels(c_terms)))) {
    cluster_part <- labels(c_terms)[grepl("cluster\\(", labels(c_terms))]

    # Extract the variable inside cluster()
    cluster_var_c <- gsub("cluster\\((.+)\\)", "\\1", cluster_part)

  }

  if (any(grepl("cluster\\(", labels(s_terms)))) {
    cluster_part <- labels(s_terms)[grepl("cluster\\(", labels(s_terms))]

    # Extract the variable inside cluster()
    cluster_var_s <- gsub("cluster\\((.+)\\)", "\\1", cluster_part)

  }

  clus <- gsub("cluster\\((.+)\\)", "\\1", cluster_part)


  # Check if both formulas have cluster terms and if they are the same
  if (!is.null(cluster_var_s) && !is.null(cluster_var_c) && cluster_var_s != cluster_var_c) {
    stop(paste("Mismatch in clustering variables: s_formula uses", cluster_var_s,
               "while c_formula uses", cluster_var_c, ". They must be the same!"))
  }


  # Check if trt is binary (0/1)
  if (!all(data[[trt]] %in% c(0, 1))) {
    stop(paste("The treatment variable is not binary"))
  }



  # Check if method is either "marginal" or "frailty"
  if (!method %in% c("marginal", "frailty")) {
    stop('Error: method must be either "marginal" or "frailty".')
  }

  if (is.null(res_time)) stop("Error: res_time cannot be NULL.")

  # Check if time is within the maximum time in data
  max_time <- max(data[[s_lhs[1]]], na.rm = TRUE)
  if (any(res_time > max_time)) stop(paste("Error: Some values in time exceed max time (", max_time, ") in data."))


  # First, check if c_lhs[2] and s_lhs[2] refer to the same variable
  if (s_lhs[2] != c_lhs[2]) {
    stop(paste("Error: The event variable in s_lhs and c_lhs must be the same.",
               "s_lhs[2] =", s_lhs[2], "but c_lhs[2] =", c_lhs[2]))
  }

  # Extract the event variable (since we confirmed they are the same)
  event_var <- s_lhs[2]

  # Extract unique values and sort them
  unique_vals <- unique(data[[event_var]])
  expected_vals <- 0:max(unique_vals)

  # Check if all required values exist
  if (!all(expected_vals %in% unique_vals)) {
    stop(paste("Error: The event variable", event_var,
               "must include all integers from 0 to", max(unique_vals),
               "without missing any. Found:", paste(unique_vals, collapse = ", ")))
  }


  ## calculate estimator

  if(all(unique_vals %in% c(0, 1))){
    # Check if est is either "surv" or "rmst"
    if (!est %in% c("surv", "rmst")) {
      stop('Error: est must be either "surv" or "rmst".')
    }

    # Check if scale is either "RD" or "RR"
    if (!scale %in% c("RD", "RR")) {
      stop('Error: scale must be either "RD" or "RR".')
    }

    est <- drs_p_est(s_formula, c_formula, trt, probs, estimand, scale,
                     method, res_time, data, id)
    var <- dr_pvar(s_formula, c_formula, trt, probs,
                   estimand, scale, method,state="s", res_time,data,id)

  }else{

    est <- drm_p_est(s_formula,c_formula,trt,probs,method,res_time,data,id)

    var <- dr_pvar(s_formula, c_formula, trt, probs,
                   estimand, scale, method,state="m", res_time,data,id)

    colnames(est$s_c) <- colnames(est$s_i) <- c("1","0","effect")

  }

  estimand <- ifelse(all(unique_vals %in% c(0, 1)), estimand,"rmt-if")
  scale <- ifelse(all(unique_vals %in% c(0, 1)),scale ,"rmt-if")

  obj <- list(estimator=est, variance=var,estimand=estimand,scale=scale,time=res_time,
              M=clus)
  class(obj) <- "dr_surv"

  return(obj)




}

#' Summary for Doubly Robust Survival Estimates
#'
#' This function provides a summary of survival estimates from `drsurv_crt()`,
#' including confidence intervals based on the specified significance level.
#'
#' @param object A list returned from `drsurv_crt()`, containing survival estimates.
#' @param alpha A numeric value (default `0.05`) specifying the significance level
#'   for confidence intervals (CI). The function calculates `(1 - alpha)%` CI.
#'
#' @return A list with summary statistics:
#'   - `$summary_s_c`: Summary of **cluster-level survival estimates** with confidence intervals.
#'   - `$summary_s_i`: Summary of **individual-level survival estimates** with confidence intervals.
#'
#' @export
summary.dr_surv <- function(object,alpha=0.05) {
  cat("Summary of doubly robust estimators:\n")
  cat("----------------------\n")
  cat("Estimand is " ,object$estimand,"with scale: ", object$scale, "\n")
  cat("The cluster-level estimators are: \n")
  tb_c <- data.frame(time=object$time,estimand = object$estimator$s_c[,3], se = sqrt(object$variance$s_c_var[,3]) )
  tb_i <- data.frame(time=object$time,estimand = object$estimator$s_i[,3], se = sqrt(object$variance$s_i_var[,3]) )

  tb_c$lower <- tb_c$estimand - qnorm(1-alpha/2)*tb_c$se
  tb_c$upper <- tb_c$estimand + qnorm(1-alpha/2)*tb_c$se

  tb_i$lower <- tb_i$estimand - qnorm(1-alpha/2)*tb_i$se
  tb_c$upper <- tb_i$estimand + qnorm(1-alpha/2)*tb_i$se

  rownames(tb_c) <- rownames(tb_i) <- NULL
  print(tb_c)
  cat("The individual-level estimators are: \n")
  print(tb_c)

  return(list(tb_i,tb_c))


}
