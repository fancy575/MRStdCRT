#' Internal Function: Doubly Robust Survival Estimation
#'
#' Estimates survival curves using marginal or frailty methods.
#' This function is for internal use only.
#'
#' @noRd

DR_S_est <- function(data,s_formula, c_formula, probs,
                     trt,times,method="marginal"){
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




  s_cov <- setdiff(all.vars(s_formula), c(rel_s[1:2],clus))
  c_cov <- setdiff(all.vars(c_formula), c(rel_c[1:2],clus))


  data <- data[order(data[,s_lhs[1]]),]
  times <- times[order(times)]


  sub_dat1 <- data[which(data[,trt]==1),]
  sub_dat0 <- data[which(data[,trt]==0),]

  if(method == "marginal"){
    s_form <- as.formula(gsub("\\+?\\s*cluster\\(.*?\\)", "", deparse(s_formula)))
    c_form <- as.formula(gsub("\\+?\\s*cluster\\(.*?\\)", "", deparse(c_formula)))

    s_fit1 <- coxph(s_form, data = sub_dat1)
    s_fit0 <- coxph(s_form, data = sub_dat0)

    c_fit1 <- coxph(c_form, data = sub_dat1)
    c_fit0 <- coxph(c_form, data = sub_dat0)


    S_est <- DR_est_mg(ftime= data[[s_lhs[1]]],delta = data[[s_lhs[2]]],trt=data[[trt]],
                          strata = data[[clus]],
                          trt_prob1 = probs[1],trt_prob0 = probs[2],
                          censor_cov = as.matrix(data[,c_cov]),
                          surv_cov = as.matrix(data[,s_cov]),
                          censor_fit1 = c_fit1$coefficients,censor_fit0 = c_fit0$coefficients,
                          surv_fit1 = s_fit1$coefficients, surv_fit0 = s_fit0$coefficients,
                          e_time = times,RMST_cal = F)


  }else if(method =="frailty"){

    s_form <- s_formula
    c_form <- c_formula

    s_fit1 <- emfrail(s_form,data=sub_dat1,distribution = emfrail_dist(dist = "gamma", theta = 2))
    s_fit0 <- emfrail(s_form,data=sub_dat0,distribution = emfrail_dist(dist = "gamma", theta = 2))

    c_fit1 <- emfrail(c_form,data=sub_dat1,distribution = emfrail_dist(dist = "gamma", theta = 2))
    c_fit0 <- emfrail(c_form,data=sub_dat0,distribution = emfrail_dist(dist = "gamma", theta = 2))


    S_est <- DR_est_fy(ftime= data[[s_lhs[1]]],delta = data[[s_lhs[2]]],trt=data[[trt]],
                          strata = data[[clus]],
                          trt_prob1 = probs[1],trt_prob0 = probs[2],
                          censor_cov = as.matrix(data[,c_cov]),
                          surv_cov = as.matrix(data[,s_cov]),
                          censor_fit1 = c_fit1$coefficients,censor_fit0 = c_fit0$coefficients,
                          beta_c1 = exp(c_fit1$logtheta) , beta_c0 = exp(c_fit0$logtheta),
                          surv_fit1 = s_fit1$coefficients, surv_fit0 = s_fit0$coefficients,
                          beta_s1 = exp(s_fit1$logtheta) ,  beta_s0 = exp(s_fit0$logtheta),
                          e_time = times,RMST_cal = F)
  }


  return(S_est)




}


#' Monotonize Survival Function (Internal)
#'
#' This function ensures a survival function is strictly decreasing.
#'
#' @noRd

S_monotone <- function(surv_vec) {
  if (!is.numeric(surv_vec)) {
    stop("Input must be a numeric vector.")
  }
  return(cummin(surv_vec))
}


#' Trapezoidal Integration for Survival Functions (Internal)
#'
#' Computes the trapezoidal integral given a series of survival values and time points.
#'
#' @noRd

trapz_integral <- function(S, time, tau) {
  if (length(S) != length(time)) stop("S and time must have the same length.")
  if (!is.numeric(S) || !is.numeric(time) || !is.numeric(tau)) stop("S, time, and tau must be numeric.")
  if (tau < min(time)) stop("Tau must be at least the minimum time.")
  if (tau > max(time)) stop("Tau exceeds the last time point in the dataset.")

  # Ensure the first time point is 0, and the first S value is 1
  if (time[1] != 0) {
    time <- c(0, time)
    S <- c(1, S)
  }


  # Interpolate S(tau) if needed
  if (!(tau %in% time)) {
    S_tau <- approx(time, S, xout = tau, rule = 2)$y
    time <- c(time[time < tau], tau)
    S <- c(S[which(time < tau)], S_tau)
  }

  # Compute trapezoidal integral
  return(pracma::trapz(time, S))
}









