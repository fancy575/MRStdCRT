#' Survival Probability Estimation using Doubly Robust Methods for single state
#'
#' This function estimates survival probabilities and restricted mean survival time using marginal or frailty working models.
#'
#' @param s_formula A survival formula for the outcome model (e.g., `Surv(time, event) ~ covariates`).
#' @param c_formula A survival formula for the censoring model.
#' @param trt A string specifying the treatment variable name.
#' @param probs A numeric vector of treatment probabilities (e.g., `c(0.5, 0.5)`).
#' @param est A string specifying the estimation method, either `"surv"` (survival) or `"rmst"` (restricted mean survival time).
#' @param scale A string specifying the scale, either `"RD"` (risk difference) or `"RR"` (risk ratio).
#' @param method A string specifying the model method, either `"marginal"` or `"frailty"`.
#' @param res_time A numeric vector specifying time points at which to estimate survival probabilities.
#' @param data A data frame containing the survival data.
#' @param id A string specifying the clustering/grouping variable name (if applicable).
#'
#' @return A list with two data frames: `s_c` (cluster-level estimates) and `s_i` (individual-level estimates).
#'
#' @import dplyr
#' @import pracma
#' @import survival
#' @import frailtyEM
#' @import abind
#' @importFrom survival
#' @export

drs_p_est <- function(s_formula, c_formula, trt, probs,
                       est, scale="RD", method, res_time,data,id){


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

  # Check if event is binary (0/1)
  if (!all(data[[s_lhs[2]]] %in% c(0, 1))) {
    stop(paste("The event variable", s_lhs[2], "must be binary (0/1)."))
  }
  if (!all(data[[c_lhs[2]]] %in% c(0, 1))) {
    stop(paste("The event variable", c_lhs[2], "must be binary (0/1)."))
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


  # Check if both formulas have cluster terms and if they are the same
  if (!is.null(cluster_var_s) && !is.null(cluster_var_c) && cluster_var_s != cluster_var_c) {
    stop(paste("Mismatch in clustering variables: s_formula uses", cluster_var_s,
               "while c_formula uses", cluster_var_c, ". They must be the same!"))
  }

  clus <- unique(c(cluster_var_s,cluster_var_c))

  # Check if trt is binary (0/1)
  if (!all(data[[trt]] %in% c(0, 1))) {
    stop(paste("The treatment variable is not binary"))
  }


  # Check if est is either "surv" or "rmst"
  if (!est %in% c("surv", "rmst")) {
    stop('Error: est must be either "surv" or "rmst".')
  }

  # Check if scale is either "RD" or "RR"
  if (!scale %in% c("RD", "RR")) {
    stop('Error: scale must be either "RD" or "RR".')
  }

  # Check if method is either "marginal" or "frailty"
  if (!method %in% c("marginal", "frailty")) {
    stop('Error: method must be either "marginal" or "frailty".')
  }

  if (is.null(res_time)) stop("Error: res_time cannot be NULL.")

  # Check if time is within the maximum time in data
  max_time <- max(data[[s_lhs[1]]], na.rm = TRUE)
  if (any(res_time > max_time)) stop(paste("Error: Some values in time exceed max time (", max_time, ") in data."))

  times <- unique(data[[s_lhs[1]]])
  ## estimate the survival function

  surv_est <- DR_S_est(data,s_formula, c_formula, probs,
                       trt,times,method)


  surv_est$S_cluster[,1:2] <- apply(surv_est$S_cluster[,1:2],2,S_monotone)
  surv_est$S_individual[,1:2] <- apply(surv_est$S_individual[,1:2],2,S_monotone)

  if(est == "rmst"){

    s_c <- vapply(res_time, function(tau) {
      apply(surv_est$S_cluster, 2, function(S_col) {
        trapz_integral(S_col, surv_est$event_time, tau)
      })
    }, numeric(ncol(surv_est$S_cluster)))


    s_i <- vapply(res_time, function(tau) {
      apply(surv_est$S_individual, 2, function(S_col) {
        trapz_integral(S_col, surv_est$event_time, tau)
      })
    }, numeric(ncol(surv_est$S_individual)))

    s_c <- t(s_c)
    s_i <- t(s_i)


  }else if(est == "surv"){

    pos <- findInterval(res_time, c(surv_est$event_time))

    s_c <- surv_est$S_cluster[pos,]

    s_i <- surv_est$S_individual[pos,]




  }


  if(scale == "RR"){
    s_c[,3] <- s_c[,1]/s_c[,2]
    s_i[,3] <- s_i[,1]/s_i[,2]
  }

  rownames(s_c) <- rownames(s_i) <- res_time
  colnames(s_c) <- colnames(s_i) <- c("S1","S0","S_diff")


  return(list(s_c = s_c, s_i = s_i))


}

#' Restricted mean time in favor (rmt-if) Estimation using Doubly Robust Methods for Multi-state
#'
#' This function estimates rmt-if estimands for cluster-level and individual-level using marginal or frailty working models.
#'
#' @param s_formula A survival formula for the outcome model (e.g., `Surv(time, event) ~ covariates`).
#' @param c_formula A survival formula for the censoring model.
#' @param trt A string specifying the treatment variable name.
#' @param probs A numeric vector of treatment probabilities (e.g., `c(0.5, 0.5)`).
#' @param method A string specifying the model method, either `"marginal"` or `"frailty"`.
#' @param res_time A numeric vector specifying time points at which to estimate survival probabilities.
#' @param data A data frame containing the survival data.
#' @param id A string specifying the clustering/grouping variable name (if applicable).
#'
#' @return A list with two data frames: `s_c` (cluster-level estimates) and `s_i` (individual-level estimates).
#'
#' @export


drm_p_est <- function(s_formula, c_formula, trt, probs,
                      method, res_time,data,id){


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


  if(s_lhs[1] != c_lhs[1]){
    stop(paste("The time variables for outocme and censoring should be the same"))

  }


  # Check if time is numeric
  if (!is.numeric(data[[s_lhs[1]]]) || !is.numeric(data[[c_lhs[1]]]) ) {
    stop(paste("The time variable must be numeric"))
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


  # Check if both formulas have cluster terms and if they are the same
  if (!is.null(cluster_var_s) && !is.null(cluster_var_c) && cluster_var_s != cluster_var_c) {
    stop(paste("Mismatch in clustering variables: s_formula uses", cluster_var_s,
               "while c_formula uses", cluster_var_c, ". They must be the same!"))
  }

  clus <- unique(c(cluster_var_s,cluster_var_c))

  # Check if trt is binary (0/1)
  if (!all(data[[trt]] %in% c(0, 1))) {
    stop(paste("The treatment variable is not binary"))
  }



  lhs_s <- as.character(s_formula[[2]])
  rhs_s <- s_formula[[3]]
  lhs_s[3] <- "event"

  lhs_mod <- as.call(c(as.name(lhs_s[1]), as.name(lhs_s[2]), as.name(lhs_s[3])))
  s_form <- as.formula(call("~", lhs_mod, rhs_s))

  lhs_c <- as.character(c_formula[[2]])
  rhs_c <- c_formula[[3]]
  lhs_c[3] <- "event"

  lhc_mod <- as.call(c(as.name(lhs_c[1]), as.name(lhs_c[2]), as.name(lhs_c[3])))
  c_form <- as.formula(call("~", lhc_mod, rhs_c))

  max_s <- max(unique_vals)
  times <- unique(data[[s_lhs[1]]])

  ## estimate the survival function
  S_c <- array(dim = c(length(times),2,max_s ))
  S_i <- array(dim = c(length(times),2,max_s ))


  for(s in 1:max_s){
    sub_dat <- data %>%
      filter(.data[[s_lhs[2]]] == 0 | .data[[s_lhs[2]]] >= s) %>%  # Filter rows where status == 0 or status >= s
      mutate(event = if_else(.data[[s_lhs[2]]] != 0, 1, 0)) %>%  # Create event column after filtering
      group_by(.data[[id]]) %>%  # Group by the ID column
      arrange(.data[[s_lhs[1]]]) %>%  # Sort by status and time
      slice(1) %>%  # Keep the first row in each group
      ungroup()  # Ungroup the data

    ## Calculate the estimate and variance for trt=1

    sub_dat <- as.data.frame(sub_dat)
    sub_dat <- sub_dat[order(sub_dat[,s_lhs[1]]),]


    surv_est <- DR_S_est(sub_dat,s_form, c_form, probs,
                         trt,times,method)

    surv_est$S_cluster[,1:2] <- apply(surv_est$S_cluster[,1:2],2,S_monotone)
    surv_est$S_individual[,1:2] <- apply(surv_est$S_individual[,1:2],2,S_monotone)



    S_c[,,s] <- surv_est$S_cluster[,1:2]
    S_i[,,s] <- surv_est$S_individual[,1:2]

  }


  # Initialize matrices to store results
  t1_matrix_c <- matrix(NA, nrow = length(res_time), ncol = max_s)
  t0_matrix_c <- matrix(NA, nrow = length(res_time), ncol = max_s)

  # Loop over each slice (q) in S_c
  for (q in 1:max_s) {
    if (q < max_s) {
      # Compute for q < max slice index
      S1 <- S_c[,1,q] * S_c[,2,q+1] - S_c[,1,q] * S_c[,2,q]
      S2 <- S_c[,2,q] * S_c[,1,q+1] - S_c[,2,q] * S_c[,1,q]
    } else {
      # Compute for the last slice (q == max)
      S1 <- S_c[,1,q] - S_c[,1,q] * S_c[,2,q]
      S2 <- S_c[,2,q] - S_c[,2,q] * S_c[,1,q]
    }

    # Apply trapz_integral for each tau in res_time
    t1_matrix_c[,q] <- vapply(res_time, function(tau) {
      trapz_integral(S1, surv_est$event_time, tau)
    }, numeric(1))

    t0_matrix_c[,q] <- vapply(res_time, function(tau) {
      trapz_integral(S2, surv_est$event_time, tau)
    }, numeric(1))
  }


  # Initialize matrices to store results
  t1_matrix_i <- matrix(NA, nrow = length(res_time), ncol = max_s)
  t0_matrix_i <- matrix(NA, nrow = length(res_time), ncol = max_s)


  # Loop over each slice (q) in S_c
  for (q in 1:max_s) {
    if (q < max_s) {
      # Compute for q < max slice index
      S1 <- S_i[,1,q] * S_i[,2,q+1] - S_i[,1,q] * S_i[,2,q]
      S2 <- S_i[,2,q] * S_i[,1,q+1] - S_i[,2,q] * S_i[,1,q]
    } else {
      # Compute for the last slice (q == max)
      S1 <- S_i[,1,q] - S_i[,1,q] * S_i[,2,q]
      S2 <- S_i[,2,q] - S_i[,2,q] * S_i[,1,q]
    }

    # Apply trapz_integral for each tau in res_time
    t1_matrix_i[,q] <- vapply(res_time, function(tau) {
      trapz_integral(S1, surv_est$event_time, tau)
    }, numeric(1))

    t0_matrix_i[,q] <- vapply(res_time, function(tau) {
      trapz_integral(S2, surv_est$event_time, tau)
    }, numeric(1))
  }





  rownames(t0_matrix_c) <-rownames(t1_matrix_c) <-  rownames(t0_matrix_i) <- rownames(t1_matrix_i) <- res_time
  colnames(t0_matrix_c) <- colnames(t1_matrix_c) <- colnames(t0_matrix_i) <- colnames(t1_matrix_i) <- 1:max_s

  s_c <- data.frame(rowSums(t1_matrix_c),rowSums(t0_matrix_c), rowSums(t1_matrix_c) -rowSums(t0_matrix_c) )
  s_i <- data.frame(rowSums(t1_matrix_i),rowSums(t0_matrix_i), rowSums(t1_matrix_i) -rowSums(t0_matrix_i) )

  rownames(s_c) <- rownames(s_i) <- res_time
  colnames(s_c) <- colnames(s_i) <- c("S1","S0","S_diff")


  return(list(s_c = s_c, s_i = s_i))


}


#' Internal Function: Jackknife resampling method for the variance calculation
#'
#'
#' @noRd

dr_pvar <- function(s_formula, c_formula, trt, probs,
                    est, scale="RD", method,state="s", res_time,data,id){

  s_terms <- terms(s_formula)

  s_lhs <- all.vars(s_formula[[2]])  # Extracts Surv(time, event)


  cluster_part <- labels(s_terms)[grepl("cluster\\(", labels(s_terms))]

    # Extract the variable inside cluster()
  clus <- gsub("cluster\\((.+)\\)", "\\1", cluster_part)



  # Check if time is within the maximum time in data
  max_time <- max(data[[s_lhs[1]]], na.rm = TRUE)
  if (any(res_time > max_time)) stop(paste("Error: Some values in time exceed max time (", max_time, ") in data."))

  times <- unique(data[[s_lhs[1]]])
  M_clus <- length(unique(data[[clus]]))
  en_clus <- unique(data[[clus]])

  jack_c <- array(dim=c(length(res_time),3,M_clus))
  jack_i <- array(dim=c(length(res_time),3,M_clus))
  if(state=="s"){
    # Check if est is either "surv" or "rmst"
    if (!est %in% c("surv", "rmst")) {
      stop('Error: est must be either "surv" or "rmst".')
    }

    # Check if scale is either "RD" or "RR"
    if (!scale %in% c("RD", "RR")) {
      stop('Error: scale must be either "RD" or "RR".')
    }
    for(m in 1:M_clus){
      lv_dt <- data[which(data[,clus] != en_clus[m] ),]
      jack_fit <- drs_p_est(s_formula, c_formula, trt, probs,
                            est, scale, method, res_time,lv_dt,id)
      jack_c[,,m] <- as.matrix(jack_fit$s_c)
      jack_i[,,m] <- as.matrix(jack_fit$s_i)

    }



  }else if (state == "m"){
    for(m in 1:M_clus){
      lv_dt <- data[which(data[,clus] != en_clus[m] ),]
      jack_fit <- drm_p_est(s_formula, c_formula, trt, probs,
                            method, res_time,lv_dt,id)
      jack_c[,,m] <- as.matrix(jack_fit$s_c)
      jack_i[,,m] <- as.matrix(jack_fit$s_i)

    }

  }


  mean_jack_c <- array(apply(jack_c, c(1,2), mean),dim=dim(jack_c))
  jack_var_c <- (M_clus / (M_clus - 1)) * apply((jack_c - mean_jack_c)^2, c(1,2), sum)

  mean_jack_i <- array(apply(jack_i, c(1,2), mean),dim=dim(jack_i))
  jack_var_i <- (M_clus / (M_clus - 1)) * apply((jack_i - mean_jack_i)^2, c(1,2), sum)

  rownames(jack_var_c) <- rownames(jack_var_i) <- res_time
  colnames(jack_var_c) <- colnames(jack_var_i) <- c("S1","S0","S_diff")

  return(list(s_c_var = jack_var_c,s_i_var = jack_var_i))

}




