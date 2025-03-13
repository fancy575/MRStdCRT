#' Model-robust Standardization Estimators for the Cluster Randomized Trials
#'
#' This function performs cluster randomized trials (CRT) analysis using model-robust standardization estimators to estimate the cluster-average and individual-average treatment effect.
#' It handles different outcome mean models (LM, LMM, GEE, GLMM) and supports both continuous, binary, and count outcomes with options for different correlation structures and scales (risk different, risk ratio and odds ratio).
#'
#' @param formula A formula for the outcome mean model, including covariates.
#' @param data A data frame where categorical variables should already be converted to dummy variables.
#' @param clus_id A string representing the column name of the cluster ID in the data frame.
#' @param trt A string representing the column name of the treatment assignment per cluster (0=control, 1=treatemnt).
#' @param prob A vector of treatment probabilities per cluster, conditional on covariates.
#' @param method A string specifying the outcome mean model. Possible values are:
#'     - 'LM': linear model on cluster-level means (continuous outcome).
#'     - 'LMM': linear mixed model on individual-level observations (continuous outcome).
#'     - 'GEE': marginal models fitted by generalized estimating equations.
#'     - 'GLMM': generalized linear mixed model.
#' @param family The link function for the outcome. Can be one of the following:
#'     - `gaussian(link = "identity")`: for continuous outcomes. Default is gaussian("identity").
#'     - `binomial(link = "logit")`: for binary outcomes.
#'     - `poisson(link = "log")`: for count outcomes.
#' @param corstr A string specifying the correlation structure for GEE models (e.g., "exchangeable", "independence").
#' @param scale A string specifying the risk measure of interest. Can be 'RD' (risk difference), 'RR' (relative risk), or 'OR' (odds ratio).
#' @param jack A numeric value (1, 2, or 3) specifying the type of jackknife standard error estimate. Type 1 is the standard jackknife, and type 3 is recommended for small numbers of clusters. Default is 1.
#' @param alpha A numeric value for the type-I error rate. Default is 0.05.
#'
#' @return A list with the following components:
#'   - `estimate`: A summary table of estimates.
#'   - `m`: Number of clusters.
#'   - `N`: Total number of observations per cluster.
#'   - `family`: The family used for the model.
#'   - `model`: The method used for the outcome mean model.
#' @import dplyr
#' @import geepack
#' @import lme4
#' @import nlme
#' @importFrom geepack geeglm
#' @importFrom stats formula
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ppact_prob <- ppact %>%
#'   group_by(CLUST) %>%
#'  mutate(first_trt = first(INTERVENTION)) %>%
#'  ungroup() %>%
#'  mutate(prob_A_1 = mean(first_trt == 1, na.rm = TRUE),  # Proportion trt=1
#'         prob_A_0 = mean(first_trt == 0, na.rm = TRUE)) %>%
#'  mutate(assigned_value = ifelse(INTERVENTION == 1, prob_A_1, prob_A_0))
#'
#'prob <- ppact_prob$assigned_value
#'
#' example <- MRStdCRT_fit(
#'   formula = PEGS ~ AGE + FEMALE + comorbid + Dep_OR_Anx + pain_count+PEGS_bl +
#'    BL_benzo_flag + BL_avg_daily + satisfied_primary + cluster(n),
#'  data = ppact,
#'  clus_id = "CLUST",
#'  trt = "INTERVENTION",
#'  prob = prob,
#'  method = "GEE",
#'  corstr = "independence",
#'  scale = "RR"
#')
#' }



MRStdCRT_fit <- function(formula, data, clus_id, trt, prob, method, family = gaussian(link="identity"),
                       corstr, scale, jack = 1, alpha=0.05){

  ################################################################
  #                                                              #
  #   Input:                                                     #
  #   formula: formula of outcome mean model.                    #
  #   data: data frame, where the categorical variables          #
  #       should be already converted to dummy variables.        #
  #   clus_id: a character of the variable name of the           #
  #            cluster id.                                       #
  #   trt: a character of the variable name of the treatment     #
  #        assignment per cluster.                               #
  #   prob: a vector of treatment probabilities per cluster      #
  #         conditional on covariates.                           #
  #   method: specifications of outcome mean models;             #
  #     potential values are:                                    #
  #       (1) 'LM' (linear model on cluster-level means),        #
  #       (2) 'LMM' (linear mixed model on individual-level      #
  #                 observations),                               #
  #       (3) 'GLMM' (marginal models fitted by generalized      #
  #                   linear mixed model),                       #
  #       (4) 'GEE' (generalized estimating equations).          #
  #   family: the link function for outcome                      #
  #           gaussian(link = "identity")                        #
  #           binomial(link = "logit")                           #
  #           poisson(link = "log").                             #
  #   corstr: correlation structure for GEE model                #
  #           "exchangeable","independence", etc.                #
  #   scale: risk differences ('RD'), relative risks ('RR'),     #
  #          and odds ratios ('OR').                             #
  #   jack: types of jackknife standard error estimates,         #
  #         potential values are 1, 2 and 3. Type 1 is           #
  #         the standard jackknife and type 3 is                 #
  #         recommended for small number of clusters.            #
  #    alpha: type-I error rate.                                 #
  ################################################################

  temp <- MRStdCRT_point(formula, data, clus_id, trt, prob,
                           family,
                           corstr, method, scale)

  data1 <- temp[[1]]
  data_clus <- temp[[2]]
  m <- nrow(data_clus)
  pes <- temp[[3]]
  cluster_names <- unique(data_clus$clus_id)
  point_est_jack <- matrix(NA, nrow = m, ncol = 3)
  for (i in cluster_names){
    data_jack <- data %>% filter(data[[clus_id]] != i)
    point_est_jack[which(cluster_names==i),] <- MRStdCRT_point(formula, data_jack, clus_id, trt, prob[which(data[[clus_id]]!=i)],
                                                                 family,
                                                                 corstr, method, scale)[[3]]
  }
  # Jackknife standard error estimation
  if (jack == 1) {
    # Implement standard jackknife
    jackse <- sqrt(apply(point_est_jack, 2, var)*(m-1)^2/m)
  } else if (jack == 2) {
    # Implement alternative jackknife method 2
    rep_point <- as.matrix(do.call(rbind, replicate(m, pes, simplify = FALSE)))
    jackse <-  as.numeric(sqrt(colSums((point_est_jack-rep_point)^2)*(m-1)/m))
  } else if (jack == 3) {
    # Implement alternative jackknife method 3
    rep_point <- as.matrix(do.call(rbind, replicate(m, pes, simplify = FALSE)))
    jackse2 <-  as.numeric(sqrt(colSums((point_est_jack-rep_point)^2)*(m-1)/m))
    jackse <- as.numeric(jackse2*sqrt(m/(m-1)))
  } else {
    stop("Invalid jackknife type specified.")
  }
  thres <- qt(p = 1-alpha/2, df = m-1)
  ciu <- pes+thres*jackse
  cil <- pes-thres*jackse
  table <- cbind(pes, jackse, cil, ciu)[-3,]
  rownames(table) <-  c("cATE",
                        "iATE")
  colnames(table) <- c("Estimate",
                       "Std, Error",
                       "CI lower",
                       "CI upper")
  table <- as.data.frame(round(table, 3))
  #test statistic for NICS
  test_sta <- pes[3]/jackse[3]
  p_val <- min((1-pt(test_sta, df = m-1, ncp = 0)),pt(test_sta, df = m-1, ncp = 0))*2
  table$`z value` <- c(round(test_sta, 3), NA)
  table$`p value` <- c(round(p_val, 3), NA)
  table[is.na(table)] <- ""

  return(list(estimate=table,m=m,N=temp[[2]]$N,
              family=family,model=method))
}
