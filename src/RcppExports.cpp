// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rate2scale
double rate2scale(double rate, double shape);
RcppExport SEXP _baclava_rate2scale(SEXP rateSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(rate2scale(rate, shape));
    return rcpp_result_gen;
END_RCPP
}
// update_scales
List update_scales(List theta);
RcppExport SEXP _baclava_update_scales(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_scales(theta));
    return rcpp_result_gen;
END_RCPP
}
// dloglik_sojourn_P_List
List dloglik_sojourn_P_List(List data_objects, List age_at_tau_hp_hats, List indolents, List theta);
RcppExport SEXP _baclava_dloglik_sojourn_P_List(SEXP data_objectsSEXP, SEXP age_at_tau_hp_hatsSEXP, SEXP indolentsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type age_at_tau_hp_hats(age_at_tau_hp_hatsSEXP);
    Rcpp::traits::input_parameter< List >::type indolents(indolentsSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(dloglik_sojourn_P_List(data_objects, age_at_tau_hp_hats, indolents, theta));
    return rcpp_result_gen;
END_RCPP
}
// dloglik_screens_List
List dloglik_screens_List(List data_objects, List age_at_tau_hp_hats, List theta);
RcppExport SEXP _baclava_dloglik_screens_List(SEXP data_objectsSEXP, SEXP age_at_tau_hp_hatsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type age_at_tau_hp_hats(age_at_tau_hp_hatsSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(dloglik_screens_List(data_objects, age_at_tau_hp_hats, theta));
    return rcpp_result_gen;
END_RCPP
}
// compute_prob_tau_obj
List compute_prob_tau_obj(List data_object, List theta, double t0);
RcppExport SEXP _baclava_compute_prob_tau_obj(SEXP data_objectSEXP, SEXP thetaSEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(compute_prob_tau_obj(data_object, theta, t0));
    return rcpp_result_gen;
END_RCPP
}
// compute_prob_tau_List
List compute_prob_tau_List(List data_objects, List theta, double t0);
RcppExport SEXP _baclava_compute_prob_tau_List(SEXP data_objectsSEXP, SEXP thetaSEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(compute_prob_tau_List(data_objects, theta, t0));
    return rcpp_result_gen;
END_RCPP
}
// rprop_age_at_tau_hp_hat_obj
NumericVector rprop_age_at_tau_hp_hat_obj(List data_object, List prob_tau, List theta, double t0);
RcppExport SEXP _baclava_rprop_age_at_tau_hp_hat_obj(SEXP data_objectSEXP, SEXP prob_tauSEXP, SEXP thetaSEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< List >::type prob_tau(prob_tauSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(rprop_age_at_tau_hp_hat_obj(data_object, prob_tau, theta, t0));
    return rcpp_result_gen;
END_RCPP
}
// rprop_age_at_tau_hp_hat_List
List rprop_age_at_tau_hp_hat_List(List data_objects, List prob_tau, List theta, double t0);
RcppExport SEXP _baclava_rprop_age_at_tau_hp_hat_List(SEXP data_objectsSEXP, SEXP prob_tauSEXP, SEXP thetaSEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type prob_tau(prob_tauSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(rprop_age_at_tau_hp_hat_List(data_objects, prob_tau, theta, t0));
    return rcpp_result_gen;
END_RCPP
}
// compute_prob_indolent_obj
NumericVector compute_prob_indolent_obj(List data_object, List theta, NumericVector age_at_tau_hp_hat);
RcppExport SEXP _baclava_compute_prob_indolent_obj(SEXP data_objectSEXP, SEXP thetaSEXP, SEXP age_at_tau_hp_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type age_at_tau_hp_hat(age_at_tau_hp_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_prob_indolent_obj(data_object, theta, age_at_tau_hp_hat));
    return rcpp_result_gen;
END_RCPP
}
// compute_prob_indolent_List
List compute_prob_indolent_List(List data_objects, List age_at_tau_hp_hats, List theta);
RcppExport SEXP _baclava_compute_prob_indolent_List(SEXP data_objectsSEXP, SEXP age_at_tau_hp_hatsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type age_at_tau_hp_hats(age_at_tau_hp_hatsSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_prob_indolent_List(data_objects, age_at_tau_hp_hats, theta));
    return rcpp_result_gen;
END_RCPP
}
// rprop_indolent_obj
IntegerVector rprop_indolent_obj(List data_object, NumericVector prob_indolent);
RcppExport SEXP _baclava_rprop_indolent_obj(SEXP data_objectSEXP, SEXP prob_indolentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob_indolent(prob_indolentSEXP);
    rcpp_result_gen = Rcpp::wrap(rprop_indolent_obj(data_object, prob_indolent));
    return rcpp_result_gen;
END_RCPP
}
// rprop_indolent_List
List rprop_indolent_List(List data_objects, List prob_indolents);
RcppExport SEXP _baclava_rprop_indolent_List(SEXP data_objectsSEXP, SEXP prob_indolentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type prob_indolents(prob_indolentsSEXP);
    rcpp_result_gen = Rcpp::wrap(rprop_indolent_List(data_objects, prob_indolents));
    return rcpp_result_gen;
END_RCPP
}
// dlog_prop_indolent_obj
NumericVector dlog_prop_indolent_obj(List data_object, NumericVector prob_indolent, IntegerVector indolent);
RcppExport SEXP _baclava_dlog_prop_indolent_obj(SEXP data_objectSEXP, SEXP prob_indolentSEXP, SEXP indolentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob_indolent(prob_indolentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indolent(indolentSEXP);
    rcpp_result_gen = Rcpp::wrap(dlog_prop_indolent_obj(data_object, prob_indolent, indolent));
    return rcpp_result_gen;
END_RCPP
}
// dlog_prop_latent_obj
NumericVector dlog_prop_latent_obj(List data_object, List prob_tau, List theta, NumericVector age_at_tau_hp_hat, NumericVector prob_indolent, IntegerVector indolent, double t0);
RcppExport SEXP _baclava_dlog_prop_latent_obj(SEXP data_objectSEXP, SEXP prob_tauSEXP, SEXP thetaSEXP, SEXP age_at_tau_hp_hatSEXP, SEXP prob_indolentSEXP, SEXP indolentSEXP, SEXP t0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< List >::type prob_tau(prob_tauSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type age_at_tau_hp_hat(age_at_tau_hp_hatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob_indolent(prob_indolentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indolent(indolentSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    rcpp_result_gen = Rcpp::wrap(dlog_prop_latent_obj(data_object, prob_tau, theta, age_at_tau_hp_hat, prob_indolent, indolent, t0));
    return rcpp_result_gen;
END_RCPP
}
// MCMC_cpp_internal
List MCMC_cpp_internal(List data_objects, List indolents, List prior, List age_at_tau_hp_hats, List theta, double epsilon_rate_H, NumericVector epsilon_rate_P, double epsilon_psi, double t0, int M, int thin, int M_thin, int n_obs, IntegerVector n_screen_positive_total, List adaptive, int verbose, int save_latent);
RcppExport SEXP _baclava_MCMC_cpp_internal(SEXP data_objectsSEXP, SEXP indolentsSEXP, SEXP priorSEXP, SEXP age_at_tau_hp_hatsSEXP, SEXP thetaSEXP, SEXP epsilon_rate_HSEXP, SEXP epsilon_rate_PSEXP, SEXP epsilon_psiSEXP, SEXP t0SEXP, SEXP MSEXP, SEXP thinSEXP, SEXP M_thinSEXP, SEXP n_obsSEXP, SEXP n_screen_positive_totalSEXP, SEXP adaptiveSEXP, SEXP verboseSEXP, SEXP save_latentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_objects(data_objectsSEXP);
    Rcpp::traits::input_parameter< List >::type indolents(indolentsSEXP);
    Rcpp::traits::input_parameter< List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< List >::type age_at_tau_hp_hats(age_at_tau_hp_hatsSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_rate_H(epsilon_rate_HSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type epsilon_rate_P(epsilon_rate_PSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_psi(epsilon_psiSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type M_thin(M_thinSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_screen_positive_total(n_screen_positive_totalSEXP);
    Rcpp::traits::input_parameter< List >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type save_latent(save_latentSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC_cpp_internal(data_objects, indolents, prior, age_at_tau_hp_hats, theta, epsilon_rate_H, epsilon_rate_P, epsilon_psi, t0, M, thin, M_thin, n_obs, n_screen_positive_total, adaptive, verbose, save_latent));
    return rcpp_result_gen;
END_RCPP
}
// model_comparison
NumericVector model_comparison(List data_object, List theta, double t0, int indolent_phase);
RcppExport SEXP _baclava_model_comparison(SEXP data_objectSEXP, SEXP thetaSEXP, SEXP t0SEXP, SEXP indolent_phaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data_object(data_objectSEXP);
    Rcpp::traits::input_parameter< List >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< int >::type indolent_phase(indolent_phaseSEXP);
    rcpp_result_gen = Rcpp::wrap(model_comparison(data_object, theta, t0, indolent_phase));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_baclava_rate2scale", (DL_FUNC) &_baclava_rate2scale, 2},
    {"_baclava_update_scales", (DL_FUNC) &_baclava_update_scales, 1},
    {"_baclava_dloglik_sojourn_P_List", (DL_FUNC) &_baclava_dloglik_sojourn_P_List, 4},
    {"_baclava_dloglik_screens_List", (DL_FUNC) &_baclava_dloglik_screens_List, 3},
    {"_baclava_compute_prob_tau_obj", (DL_FUNC) &_baclava_compute_prob_tau_obj, 3},
    {"_baclava_compute_prob_tau_List", (DL_FUNC) &_baclava_compute_prob_tau_List, 3},
    {"_baclava_rprop_age_at_tau_hp_hat_obj", (DL_FUNC) &_baclava_rprop_age_at_tau_hp_hat_obj, 4},
    {"_baclava_rprop_age_at_tau_hp_hat_List", (DL_FUNC) &_baclava_rprop_age_at_tau_hp_hat_List, 4},
    {"_baclava_compute_prob_indolent_obj", (DL_FUNC) &_baclava_compute_prob_indolent_obj, 3},
    {"_baclava_compute_prob_indolent_List", (DL_FUNC) &_baclava_compute_prob_indolent_List, 3},
    {"_baclava_rprop_indolent_obj", (DL_FUNC) &_baclava_rprop_indolent_obj, 2},
    {"_baclava_rprop_indolent_List", (DL_FUNC) &_baclava_rprop_indolent_List, 2},
    {"_baclava_dlog_prop_indolent_obj", (DL_FUNC) &_baclava_dlog_prop_indolent_obj, 3},
    {"_baclava_dlog_prop_latent_obj", (DL_FUNC) &_baclava_dlog_prop_latent_obj, 7},
    {"_baclava_MCMC_cpp_internal", (DL_FUNC) &_baclava_MCMC_cpp_internal, 17},
    {"_baclava_model_comparison", (DL_FUNC) &_baclava_model_comparison, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_baclava(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
