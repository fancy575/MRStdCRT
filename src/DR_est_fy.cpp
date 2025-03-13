#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]


Rcpp::List DR_est_fy(const arma::vec& ftime,const arma::vec& delta, const arma::vec trt,
                       const arma::vec& strata, const double& trt_prob1, const double& trt_prob0,
                       const arma::mat& censor_cov, const arma::mat& surv_cov,
                       const arma::vec& censor_fit1, const arma::vec& censor_fit0,
                       const double& beta_c1, const double& beta_c0,
                       const arma::vec& surv_fit1, const arma::vec& surv_fit0,
                       const double& beta_s1, const double& beta_s0,
                       const arma::vec& e_time,const bool& RMST_cal){

  arma::vec strata_M = unique(strata);
  int M = strata_M.n_elem;
  int N = ftime.n_elem;

  //arma::vec e_time = unique(ftime);



  int t = e_time.n_elem;

  arma::cube DR_est_cube(M,3,t);
  arma::vec cluster_size(M);


  // Construct counting process and risk process
  arma::mat e_time_expand = repmat(e_time.t(),N,1);
  arma::mat ftime_expand= repmat(ftime,1,t);
  arma::mat delta_expand = repmat(delta,1,t);
  arma::umat dN_pre = (ftime_expand == e_time_expand);
  arma::mat dN_pre_convert = arma::conv_to<arma::mat>::from(dN_pre);
  arma::mat dN = dN_pre_convert % delta_expand;
  arma::mat dN_c = dN_pre_convert % (1-delta_expand);

  arma::umat risk_pre = (ftime_expand >= e_time_expand);
  arma::mat risk = arma::conv_to<arma::mat>::from(risk_pre);

  // Calculate exponential part and baseline hazard
  // trt=1
  arma::mat s_exp1 = exp(surv_cov * surv_fit1);
  arma::mat c_exp1 = exp(censor_cov * censor_fit1);

  arma::uvec trt_row1 = find(trt==1);
  int N_sub1 = trt_row1.n_elem;


  arma::mat s_exp1_sub = s_exp1.rows(trt_row1);
  arma::mat c_exp1_sub = c_exp1.rows(trt_row1);

  arma::mat risk1_sub = risk.rows(trt_row1);
  arma::mat dN1_sub = dN.rows(trt_row1);
  arma::mat dN1_c_sub = dN_c.rows(trt_row1);



  // Calculate baseline hazard
  arma::mat s_exp1_expand = repmat(s_exp1,1,t);
  arma::mat c_exp1_expand = repmat(c_exp1,1,t);

  arma::mat s_exp1_sub_expand = repmat(s_exp1_sub,1,t);
  arma::mat c_exp1_sub_expand = repmat(c_exp1_sub,1,t);



  arma::mat S0_s_elem1 = s_exp1_sub_expand % risk1_sub;
  arma::mat S0_s_expand1 = repmat(sum(S0_s_elem1,0),N_sub1,1);
  arma::vec dlambda0_s1 = sum(dN1_sub/S0_s_expand1,0).t();
  //arma::mat dlambda0_s1_expand = repmat(dlambda0_s1.t(),N,1);
  arma::vec Lambda0_s1 = cumsum(dlambda0_s1);
  //arma::vec Surv0_s1 = exp(-Lambda0_s1);



  arma::mat S0_c_elem1 = c_exp1_sub_expand % risk1_sub;
  arma::mat S0_c_expand1 = repmat(sum(S0_c_elem1,0),N_sub1,1);
  arma::vec dlambda0_c1 = sum(dN1_c_sub/S0_c_expand1,0).t();
  //arma::mat dlambda0_c1_expand = repmat(dlambda0_c1.t(),N,1);
  arma::vec Lambda0_c1 = cumsum(dlambda0_c1);
  //arma::vec Surv0_c1 = exp(-Lambda0_c1);





  // Calculate martingale matrix

  // lambda_c
  arma::mat dlambda0_c_expand1 = repmat(dlambda0_c1.t(),N,1);
  arma::mat dlambda_c1 = dlambda0_c_expand1 % c_exp1_expand;
  arma::mat Lambda_c1 = cumsum(dlambda_c1,1);
  arma::mat laplace_c1 = beta_c1 / (beta_c1+Lambda_c1);



  arma::mat dM1_c = dN_c - laplace_c1 % dlambda0_c_expand1 % c_exp1_expand % risk ;


  // For Kc and H1



  arma::mat Lambda0_c1_expand = repmat(Lambda0_c1.t(),N,1);
  arma::mat Lambda0_s1_expand = repmat(Lambda0_s1.t(),N,1);

  arma::mat Kc1 = pow(beta_c1/(beta_c1 + Lambda0_c1_expand % c_exp1_expand ), beta_c1);
  arma::mat H11 = pow(beta_s1/(beta_s1 + Lambda0_s1_expand % s_exp1_expand ), beta_s1);




  // Integral part
  arma::mat integral_part1 = dM1_c / Kc1 / H11;


  // for trt =0
  arma::mat s_exp0 = exp(surv_cov * surv_fit0);
  arma::mat c_exp0 = exp(censor_cov * censor_fit0);

  arma::uvec trt_row0 = find(trt==0);
  int N_sub0 = trt_row0.n_elem;


  arma::mat s_exp0_sub = s_exp0.rows(trt_row0);
  arma::mat c_exp0_sub = c_exp0.rows(trt_row0);

  arma::mat risk0_sub = risk.rows(trt_row0);
  arma::mat dN0_sub = dN.rows(trt_row0);
  arma::mat dN0_c_sub = dN_c.rows(trt_row0);



  // Calculate baseline hazard
  arma::mat s_exp0_expand = repmat(s_exp0,1,t);
  arma::mat c_exp0_expand = repmat(c_exp0,1,t);

  arma::mat s_exp0_sub_expand = repmat(s_exp0_sub,1,t);
  arma::mat c_exp0_sub_expand = repmat(c_exp0_sub,1,t);



  arma::mat S0_s_elem0 = s_exp0_sub_expand % risk0_sub;
  arma::mat S0_s_expand0 = repmat(sum(S0_s_elem0,0),N_sub0,1);
  arma::vec dlambda0_s0 = sum(dN0_sub/S0_s_expand0,0).t();
  //arma::mat dlambda0_s0_expand = repmat(dlambda0_s0.t(),N,1);
  arma::vec Lambda0_s0 = cumsum(dlambda0_s0);
  //arma::vec Surv0_s0 = exp(-Lambda0_s0);



  arma::mat S0_c_elem0 = c_exp0_sub_expand % risk0_sub;
  arma::mat S0_c_expand0 = repmat(sum(S0_c_elem0,0),N_sub0,1);
  arma::vec dlambda0_c0 = sum(dN0_c_sub/S0_c_expand0,0).t();
  //arma::mat dlambda0_c0_expand = repmat(dlambda0_c0.t(),N,1);
  arma::vec Lambda0_c0 = cumsum(dlambda0_c0);
  //arma::vec Surv0_c0 = exp(-Lambda0_c0);





  // Calculate martingale matrix

  // lambda_c
  arma::mat dlambda0_c_expand0 = repmat(dlambda0_c0.t(),N,1);
  arma::mat dlambda_c0 = dlambda0_c_expand0 % c_exp0_expand;
  arma::mat Lambda_c0 = cumsum(dlambda_c0,1);
  arma::mat laplace_c0 = beta_c0 / (beta_c0 + Lambda_c0);

  arma::mat dM0_c = dN_c - laplace_c0  % dlambda0_c_expand0 % c_exp0_expand % risk;


  // For Kc and H1


  arma::mat Lambda0_c0_expand = repmat(Lambda0_c0.t(),N,1);
  arma::mat Lambda0_s0_expand = repmat(Lambda0_s0.t(),N,1);

  arma::mat Kc0 = pow(beta_c0/(beta_c0 + Lambda0_c0_expand % c_exp0_expand ), beta_c0);
  arma::mat H10 = pow(beta_s0/(beta_s0 + Lambda0_s0_expand % s_exp0_expand ), beta_s0);



  // Integral part
  arma::mat integral_part0 = dM0_c / Kc0 / H10;

  arma::mat trt1_expand = repmat(trt,1,t);
  arma::mat trt0_expand = repmat(1-trt,1,t);

  // Assemble element
  // trt=1

  arma::mat part1_1 = trt1_expand % risk / Kc1/trt_prob1;
  part1_1.elem(find_nonfinite(part1_1)).zeros();
  arma::mat part2_1 = (trt1_expand - trt_prob1)/trt_prob1 % H11;
  part2_1.elem(find_nonfinite(part2_1)).zeros();
  arma::mat part3_int1 = cumsum(integral_part1,1);
  arma::mat part3_1 = trt1_expand/trt_prob1 % H11 % part3_int1;
  part3_1.elem(find_nonfinite(part3_1)).zeros();

  arma::mat S_ci_1 = part1_1 - part2_1 + part3_1;



  //trt=0
  arma::mat part1_0 = trt0_expand % risk/Kc0/trt_prob0;
  part1_0.elem(find_nonfinite(part1_0)).zeros();
  arma::mat part2_0 = (trt0_expand - trt_prob0)/trt_prob0 % H10;
  part2_0.elem(find_nonfinite(part2_0)).zeros();
  arma::mat part3_int0 = cumsum(integral_part0,1);
  arma::mat part3_0 = trt0_expand/trt_prob0 % H10 % part3_int0;
  part3_0.elem(find_nonfinite(part3_0)).zeros();

  arma::mat S_ci_0 = part1_0 - part2_0 + part3_0;


  for(int m=0; m <M; m++){
    arma::uvec strata_id = find(strata == strata_M(m));
    int Ni = strata_id.n_elem;

    arma::mat temp_S_ci_1 = S_ci_1.rows(strata_id);
    arma::mat temp_S_ci_0  = S_ci_0.rows(strata_id);
    arma::mat temp_S_ci_diff = temp_S_ci_1 - temp_S_ci_0;

    arma::mat temp_Sc_1 = mean(temp_S_ci_1,0);
    arma::mat temp_Sc_0 = mean(temp_S_ci_0,0);
    arma::mat temp_Sc_diff = mean(temp_S_ci_diff,0);

    DR_est_cube(span(m,m),span(0,0),span(0,t-1)) = temp_Sc_1;
    DR_est_cube(span(m,m),span(1,1),span(0,t-1)) = temp_Sc_0;
    DR_est_cube(span(m,m),span(2,2),span(0,t-1)) = temp_Sc_diff;

    cluster_size(m) = Ni;


  }

  // Calculate cluster level and individual level S
  arma::mat out_S_cluster(t,3);
  arma::mat cluster_mean = mean(DR_est_cube,0);
  out_S_cluster(span(0,t-1),span(0,2)) = cluster_mean.t();
  out_S_cluster.for_each([](mat::elem_type& val){if(val > 1){val = 1;}else if(val <0){val=0;} } );



  arma::mat cluster_size_prep = repmat(cluster_size,1,3);
  arma::cube out_ind_prep = DR_est_cube.each_slice() % cluster_size_prep;
  arma::mat out_ind_sum = sum(out_ind_prep,0);
  int all_N = sum(cluster_size);
  arma::mat out_S_ind = out_ind_sum.t() / all_N;
  out_S_ind.for_each([](mat::elem_type& val){if(val > 1){val = 1;}else if(val <0){val=0;} } );

  // Calculate RMST


  arma::vec time_diff1 = e_time;
  arma::vec time_diff2(t);
  time_diff2(span(1,t-1)) = e_time(span(0,t-2));
  arma::vec time_diff = time_diff1 - time_diff2;
  arma::mat S_cluster_prep1 = out_S_cluster(span(0,t-1),span(0,2));
  arma::mat S_cluster_prep2 = out_S_cluster(span(0,t-2),span(0,2));
  arma::vec insert_1 = {1,1,0};
  S_cluster_prep2.insert_rows(0,insert_1.t());

  //S_cluster_prep1.for_each([](mat::elem_type& val){if(val > 1){val = 1;}else if(val <0){val=0;} } );
  //S_cluster_prep2.for_each([](mat::elem_type& val){if(val > 1){val = 1;}else if(val <0){val=0;} } );

  if(RMST_cal){
    arma::mat time_diff_expand = repmat(time_diff,1,3);
    arma::mat RMST_cluster = (S_cluster_prep2 + S_cluster_prep1) % time_diff_expand/2;


    arma::mat S_ind_prep1 = out_S_ind(span(0,t-1),span(0,2));
    arma::mat S_ind_prep2 = out_S_ind(span(0,t-2),span(0,2));
    S_ind_prep2 .insert_rows(0,insert_1.t());

    arma::mat RMST_ind = (S_ind_prep2 + S_ind_prep1) % time_diff_expand/2;


    arma::mat RMST_cluster_out = cumsum(RMST_cluster,0);
    arma::mat RMST_ind_out = cumsum(RMST_ind,0);

    return Rcpp::List::create(Rcpp::Named("S_individual") = out_S_ind,
                      Rcpp::Named("S_cluster") = out_S_cluster,
                      Rcpp::Named("Cluster_size") = cluster_size,
                      Rcpp::Named("event_time") = e_time,
                      Rcpp::Named("RMST_cluster_out") = RMST_cluster_out,
                      Rcpp::Named("RMST_ind_out") = RMST_ind_out




                            );



  }else {
    return Rcpp::List::create(Rcpp::Named("S_individual") = out_S_ind,
                      Rcpp::Named("S_cluster") = out_S_cluster,
                      Rcpp::Named("Cluster_size") = cluster_size,
                      Rcpp::Named("event_time") = e_time

                            );



  }




}
