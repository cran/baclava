#include <Rcpp.h>
#include <Rmath.h>
#include <string>
#include <RcppNumerical.h>

using namespace Numer;
using namespace Rcpp;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

// -----------------------------------------------------------------------------
// Within each parameter group, there are 3 outcome groups:
//   Grp1 satisfy I(T_i < age_at_tau_hp_i < age_at_tau_pc_i)
//   Grp2 satisfy I(age_at_tau_hp_i < T_i < age_at_tau_pc_i)
//   Grp3 satisfy I(age_at_tau_hp_i < age_at_tau_pc_i == T_i)
// I use this notation in comments.
//   Grp1 come from censored cases; 
//   Grp2 come from screen or censored cases;
//   Grp3 come from clinical cases
// -----------------------------------------------------------------------------

// Re-creation of R's rep when two vectors are provided
//
// @param x NumericVector The values to be repeated.
// @param times IntegerVector The number of times to repeat each value.
//   Vector must be of the same length as input `x`.
//
// @returns NumericVector of length sum(times)
//
// @keywords internal
NumericVector repVec(const NumericVector& x, const IntegerVector& times) {
  R_xlen_t n = times.length();
  if (n != 1 && n != x.length()) stop("Invalid 'times' value");
  R_xlen_t n_out = std::accumulate(times.begin(), times.end(), 0);
  NumericVector res = no_init(n_out);
  auto begin = res.begin();
  for (R_xlen_t i = 0, ind = 0; i < n; ind += times[i], ++i) {
      auto start = begin + ind;
      auto end = start + times[i];
      std::fill(start, end, x[i]);
  }
  return res;
}


//
// Weibull ////////


// Given rate and shape values, calculate the corresponding scale value
//
// @param rate double The rate parameter of a Weibull distribution. Must be > 0.
// @param shape double The shape parameter of a Weibull distribution. Must be > 0.
//
// @returns double The scale parameter of the Weibull distribution. Will be > 0.
//
// exported auxiliary function
// [[Rcpp::export]]
double rate2scale(double rate, double shape) {
  NumericVector temp(1, rate);
  return Rcpp::pow(temp, -1.0 / shape)[0];
}

// Given a parameter List with updated rates, calculate new scales
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The input `theta` with updated elements `$scale_H`
//   and `$scale_P`.
//
// exported auxiliary function
// [[Rcpp::export]]
List update_scales(List theta) {
  
  // all participants have the same healthy distribution
  theta["scale_H"] = rate2scale(theta["rate_H"], theta["shape_H"]);
  
  // the preclinical distribution can vary across participants
  NumericVector rate_P = theta["rate_P"];
  NumericVector scale_P = no_init(rate_P.size());
  std::transform(rate_P.begin(),
                 rate_P.end(),
                 scale_P.begin(), [=] (double x) {return rate2scale(x, theta["shape_P"]); });
  theta["scale_P"] = scale_P;
  
  return theta;
}

// Probability P(a < X <= b) for a Weibull distribution
//
// Note: Differences are taken between element i of input `b` and element i of
//   input `a`.
//
// @param a NumericVector A vector of one or more quantiles.
// @param b NumericVector A vector of one or more quantiles. Must be the same
//   length as input `a`.
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
//
// @returns NumericVector The ith element is the difference in the Weibull
//   CDF at the ith quantiles of inputs `a` and `b`.
//
// @keywords internal
NumericVector pweibull_ab(NumericVector a, 
                          NumericVector b, 
                          double shape, 
                          double scale) {
  return pweibull(b, shape, scale, true, false) -
    pweibull(a, shape, scale, true, false);
}

// A random sample from the truncated Weibull distribution.
//
// Note: Upper and lower limits are taken in pairs from inputs `a` and `b`; 
//   i.e., element i of `a` is the lower bound, element i of `b` is the upper
//   bound.
//
// @param a NumericVector A vector of one or more quantiles.
// @param b NumericVector A vector of one or more quantiles. Must be the same
//   length as input `a`.
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
//
// @returns NumericVector The generated samples.
//
// @keywords internal
NumericVector rweibull_trunc(NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale) {
    
  NumericVector low = pweibull(a, shape, scale);
  NumericVector high = pweibull(b, shape, scale);
    
  NumericVector ru = runif(low.size());
  NumericVector u = low + (high - low) * ru;
    
  return qweibull(u, shape, scale);
}

// Truncated Weibull Probability Density
//
// @param x NumericVector A vector of one or more quantiles at which the
//   density is evaluated.
// @param a NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`. The lower boundary of the reduced support.
// @param b NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`. The upper boundary of the reduced support
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
// @param uselog bool If TRUE, probabilities p are returned as log(p).
//
// @returns NumericVector If uselog = false
//     f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)}
//   if uselog = true
//     ln[ f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)} ]
//
// @keywords internal
NumericVector dweibull_trunc(NumericVector x, 
                             NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale, 
                             bool uselog) {
  
  if (uselog) {
    // ln[ f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)} ]
      return dweibull(x, shape, scale, true) - log(pweibull_ab(a, b, shape, scale));
  } else {
    // f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)}
      return dweibull(x, shape, scale, false) / pweibull_ab(a, b, shape, scale);
  }
}

//
// theta ////////

// Store a new rate_H value and update its corresponding scale_H.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param rate_H double The new rate value to be stored.
//
// @returns List The input `theta` with updated elements `$rate_H` and
//  `$scale_H`.
//
// @keywords internal
List add_rate_H(List theta, double rate_H) {
  theta["rate_H"] = rate_H;
  return update_scales(theta);
}

// Store a new rate_P vector and update the corresponding scale_P.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param rate_P NumericVector The new rate values to be stored.
//
// @returns List The input `theta` with updated elements `$rate_P` and
//  `$scale_P`.
//
// @keywords internal
List add_rate_P(List theta, NumericVector rate_P) {
  theta["rate_P"] = rate_P;
  return update_scales(theta);
}

// Store a new beta vector.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param beta NumericVector The new value to be stored in element `$beta`.
//
// @returns List The input `theta` with updated element `$beta`.
//
// @keywords internal
List add_beta(List theta, NumericVector beta) {
  theta["beta"] = beta;
  return theta;
}

// Store a new psi value.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param psi double The new value to be stored in element `$psi`.
//
// @returns List The input `theta` with updated element `$psi`.
//
// @keywords internal
List add_psi(List theta, double psi) {
  theta["psi"] = psi;
  return theta;
}

//
// Gibbs theta ////////

// Identify the number of screens that occurred after the current estimated
//   transition to preclinical.
//
// @param data_object List The data for a single censor type.
// @param age_at_tau_hp_hats NumericVector the current tau_hp estimates.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The input `theta` with updated element `$beta`.
//
// @keywords internal
IntegerVector gibbs_beta_obj(List data_object, 
                             NumericVector age_at_tau_hp_hat,
                             List theta) {
  
  NumericVector betas = theta["beta"];
  int n_screen_types = betas.size();
  
  // count the number of screens that occurred after the transition time
  IntegerVector missed(n_screen_types, 0);
  
  // count the total number of screens that occurred after the current
  // estimated age at time of healthy -> preclinical transition
  
  // extract vectorized age data for current censor type
  List age_screen = data_object["ages_screen"];
  
  // no screens means that beta is 0
  if (age_screen.size() == 0) {
    return missed;
  }
  
  NumericVector ages = age_screen["values"];
  IntegerVector screen_types = age_screen["types"];
  IntegerVector age_starts = age_screen["starts"];
  IntegerVector age_ends = age_screen["ends"];

  for (int i = 0; i < age_at_tau_hp_hat.size(); ++i) {
    for (int j = age_starts[i]; j <= age_ends[i]; ++j) {
      if (ages[j] > age_at_tau_hp_hat[i]) missed[screen_types[j]] += 1;
    }
  }
  return missed;
}

List gibbs_beta_List(List data_objects, 
                     List prior, 
                     List age_at_tau_hp_hats,
                     List theta, 
                     IntegerVector n_screen_positive) {
  
  NumericVector beta = theta["beta"];

  IntegerVector n_screen_since_transition(beta.size(), 0);
  for (int i = 0; i < data_objects.size(); ++i) {
    n_screen_since_transition += gibbs_beta_obj(data_objects[i], 
                                                age_at_tau_hp_hats[i], theta);        
  }

  NumericVector prior_a_beta = prior["a_beta"];
  NumericVector prior_b_beta = prior["b_beta"];
  NumericVector beta_new = no_init(beta.size());
  double a_beta, b_beta, a_n, b_n;
  
  for (int i = 0; i < beta.size(); ++i) {
    a_beta = prior_a_beta[i];
    b_beta = prior_b_beta[i];
    if (a_beta < 1e-12 && b_beta < 1e-12) {
      beta_new[i] = beta[i];
    } else {
      // a* = a0 + {# of screens that successfully detected preclinical cancers}
      a_n = a_beta + n_screen_positive[i];
      // b* = b0 + {# of screens that failed to detect preclinical cancers}
      b_n = b_beta + n_screen_since_transition[i] - n_screen_positive[i];
      beta_new[i] = rbeta(1, a_n, b_n)[0];
    }
  }

  return add_beta(theta, beta_new);
}

//
// Likelihood ////////

// The log-likelihood terms involving the sojourn time in the healthy state
//   for all individuals of a specific censor type
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type 
//   under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector  Each element provides
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
NumericVector dloglik_sojourn_H_obj(List data_object, 
                                    List theta,
                                    NumericVector age_at_tau_hp_hat,
                                    double t0) {
    
  NumericVector endpoint_time = data_object["endpoint_time"];
  int n = data_object["n"];
    
  NumericVector result = no_init(n);
    
  LogicalVector isInfinite = is_infinite(age_at_tau_hp_hat);
    
  // for participants that have an infinite age at time of transition, that is
  // so-called "healthy at censoring time",
  // F_h(T_i - t0, shape, scale); Pr(tau_hp <= T_i - t0)
  // S_H(c_i - t_0; \lambda_H, \alpha_H) = 1 - F_h
  NumericVector vec1 = endpoint_time[isInfinite];
  NumericVector result_infinite = pweibull(vec1 - t0,
                                           theta["shape_H"],
                                           theta["scale_H"],
                                           false, true);
  result[isInfinite] = result_infinite;
    
  // for participants that have a finite estimated age at time of transition,
  // i.e., those participants in Grp2 and Grp3
  // f_h(age_at_tau_hp_hat_i - t0, shape, scale); Pr(tau_hp == tau_hp_hat)
  vec1 = age_at_tau_hp_hat[!isInfinite];
  NumericVector result_finite = dweibull(vec1 - t0,
                                         theta["shape_H"],
                                         theta["scale_H"],
                                         true);
    
  result[!isInfinite] = result_finite;
  return result;
}

// Sum of the log-likelihood terms involving the sojourn time in the
//   healthy state over all individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type and distribution groups.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns double 
// sum_i = 1, n
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
double dloglik_sojourn_H_sum(List data_objects, 
                             List age_at_tau_hp_hats,
                             List theta, 
                             double t0) {
  
  double result = 0.0;
  for (int i = 0; i < data_objects.size(); ++i) {
    result += sum(dloglik_sojourn_H_obj(data_objects[i], theta, 
                                        age_at_tau_hp_hats[i], t0));        
  }
  
  return result;
}

// The log-likelihood terms involving the sojourn time in the progressive
//   preclinical state for all individuals of a specific censor type and
//   parameter set
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector 0/1 indicating if is not/is indolent.
//
// @returns NumericVector(data_object["n"])  Each element providing
//   Grp1 {0}
//   Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
//   Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
NumericVector dloglik_sojourn_P_obj(List data_object, 
                                    List theta, 
                                    NumericVector age_at_tau_hp_hat,
                                    IntegerVector indolent) {
    
  NumericVector endpoint_time = data_object["endpoint_time"];
    
  int endpoint_type = data_object["endpoint_type"];
  int n = data_object["n"];
    
  NumericVector result(n, 0.0);
    
  // subset includes participants with progressive tumors from Grp2 and Grp3
  LogicalVector subset = is_finite(age_at_tau_hp_hat) & (indolent == 0);
  NumericVector tau_p_hat = endpoint_time - age_at_tau_hp_hat;
    
  int irateP = data_object["irateP"];
  NumericVector scale_P = theta["scale_P"];
    
  if (endpoint_type == 3) {
    // Clinical cases
    // Grp3
    // ln[ f_p(T - age_at_tau_hp_hat)]; ln[ Pr(tau_p == tau_p_hat) ]
    result = dweibull(tau_p_hat, theta["shape_P"], scale_P[irateP], true);
  } else {
    // Censored and preclinical cases
    // Grp2
    // ln[ F_P(T - age_at_tau_hp_hat) ]; ln[ Pr(tau_p <= tau_p_hat) ]
    result = pweibull(tau_p_hat, theta["shape_P"], scale_P[irateP], false, true);
  }
  result[!subset] = 0.0;
    
  return result;
}

// [[Rcpp::export]]
List dloglik_sojourn_P_List(List data_objects, List age_at_tau_hp_hats, List indolents, List theta) {
  
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dloglik_sojourn_P_obj(data_objects[i], theta, age_at_tau_hp_hats[i],
                                      indolents[i]);
  }
  return result;  
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of being indolent for a single censor type
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param indolent IntegerVector 0/1 indicating is not/is indolent.
//
// @returns NumericVector Each element provides
//   ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
//
// @keywords internal
NumericVector dloglik_indolent_obj(List theta, IntegerVector indolent) {
  double psi = theta["psi"];
    
  int n = indolent.size();
  NumericVector result = no_init(n);
  result[indolent == 1] = log(psi);
  result[indolent == 0] = log(1.0 - psi);
  return result;    
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of being indolent for all individuals
//
// @param indolents List The current estimates for 0/1 indolent.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List 
//   Each element is a NumericVector providing
//     ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
//   for a unique endpoint_type/beta/rate_P combination
//
// @keywords internal
List dloglik_indolent_List(List indolents, List theta) {
  List result(indolents.size());
  for (int i = 0; i < indolents.size(); ++i) {
    result[i] = dloglik_indolent_obj(theta, indolents[i]);
  }
  return result;
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of a screen successfully detecting 
//   preclinical status for a single censor type
//
// @param data_object List The data for a single censor type and parameter
//   group.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type 
//   under consideration.
//
// @returns NumericVector Each element provides 
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
NumericVector dloglik_screens_obj(List data_object, 
                                  List theta, 
                                  NumericVector age_at_tau_hp_hat) {
  
  // extracting vectorized age at screening data
  List ages_screen = data_object["ages_screen"];
  
  if (ages_screen.size() == 0) {
    int n = data_object["n"];
    NumericVector tmp(n, 0.0);
    return tmp;
  }
  
  NumericVector age_screen = ages_screen["values"];
  IntegerVector age_starts = ages_screen["starts"];
  IntegerVector age_ends = ages_screen["ends"];
  IntegerVector screen_types = ages_screen["types"];
  
  // =1 if screen detected preclinical cancer
  IntegerVector n_screen_positive = data_object["n_screen_positive"];
  
  NumericVector betas = theta["beta"];
  
  int endpoint_type = data_object["endpoint_type"];

  NumericVector results(age_starts.size(), 0.0);
  int j;
  
  for (int i = 0; i < age_starts.size(); ++i) {
    // for each participant, determine if screens missed/caught a preclinical case
    j = age_starts[i];
    while (j < age_ends[i]) {
      if (age_screen[j] > age_at_tau_hp_hat[i]) {
        results[i] += std::log(1.0 - betas[screen_types[j]]);
      }
      j++;
    }
    if (endpoint_type == 1) {
      // if preclinical, last screen caught disease
      results[i] += std::log(betas[screen_types[age_ends[i]]]);
    } else {
      if (age_screen[age_ends[i]] > age_at_tau_hp_hat[i]) {
        results[i] += std::log(1.0 - betas[screen_types[age_ends[i]]]);
      }
    }
  }
  return results;
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of a screen successfully detecting 
//   preclinical status for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List A list. Each element is a NumericVector providing
//     n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//   for a unique endpoint_type/beta/rate_P combination
//
// @keywords internal
// [[Rcpp::export]]
List dloglik_screens_List(List data_objects, List age_at_tau_hp_hats, List theta) {
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dloglik_screens_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
  }
  return result;  
}

// Not used but might come in handy later
NumericVector dlog_likelihood_obj(List data_object, 
                                  List theta,
                                  NumericVector age_at_tau_hp_hat,
                                  IntegerVector indolent, 
                                  double t0) {
  
  // vector of length data_object["n"]
  NumericVector dlog_H = dloglik_sojourn_H_obj(data_object,
                                               theta,
                                               age_at_tau_hp_hat,
                                               t0);
  // vector of length data_object["n"]
  NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, 
                                               theta, 
                                               age_at_tau_hp_hat,
                                               indolent);
  
  // vector of length data_object["n"]
  NumericVector dlog_S = dloglik_screens_obj(data_object, 
                                             theta, 
                                             age_at_tau_hp_hat);

  NumericVector dlog_I = dloglik_indolent_obj(theta, indolent);
  
  return dlog_H + dlog_P + dlog_S + dlog_I;
  
}


List dlog_likelihood_List(List data_objects, 
                          List age_at_tau_hp_hats, 
                          List indolents, 
                          List theta, 
                          double t0) {
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dlog_likelihood_obj(data_objects[i], 
                                    theta,
                                    age_at_tau_hp_hats[i],
                                    indolents[i], 
                                    t0);
  }
  return result;  
}

// Not used but might come in handy later
double dlog_likelihood(List data_objects, 
                       List indolents, 
                       List age_at_tau_hp_hats, 
                       List theta,
                       double t0) {

  double result = 0.0;
  for (int i = 0; i < data_objects.size(); ++i) {
    result += sum(dlog_likelihood_obj(data_objects[i], 
                                      theta, 
                                      age_at_tau_hp_hats[i],
                                      indolents[i],
                                      t0));
  }
  
  return result;  
}



// The log-likelihood terms involving the sojourn time in the progressive 
//   preclinical state and the indolence parameters for all individuals in
//   of single censor type.
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector =1 if indolent.
//
// @returns NumericVector(data_object["n"]) Each element providing
//   Grp1 ln[psi]
//   Grp2 {I(ind_i == 1) ln[psi] + 
//         I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//   Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//   for a single censor type/parameter group combination.
// @keywords internal
NumericVector dloglik_PI_obj(List data_object, 
                             List theta, 
                             NumericVector age_at_tau_hp_hat,
                             IntegerVector indolent) {
  // Grp1 {0} +
  // ln[ F_P(tau_hp_hat_i) ] I(ind_i == 0) Grp2 +
  // ln[ f_P(tau_pc_hat_i) ] I(ind_i == 0) Grp3
  // NumericVector of length data_object["n"]
  NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta,
                                               age_at_tau_hp_hat, indolent);
    
  // ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
  // NumericVector of length data_object["n"]
  NumericVector dlog_I = dloglik_indolent_obj(theta, indolent);
    
  return dlog_P + dlog_I;
}

// Sum of the log-likelihood terms involving the sojourn time in the 
//   progressive preclinical state and the indolence parameters over all 
//   individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns NumericVector. Each element provides 
//     sum i = 1, n
//       Grp1 ln[psi]
//       Grp2 {I(ind_i == 1) ln[psi] + 
//             I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//       Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//   for a unique endpoint type/beta/rate_P combination
// @keywords internal
NumericVector dloglik_PI_sum(List data_objects, 
                             List indolents, 
                             List age_at_tau_hp_hats, 
                             List theta) {
  
  NumericVector result(data_objects.size(), 0.0);
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] += sum(dloglik_PI_obj(data_objects[i], theta, 
                                    age_at_tau_hp_hats[i], indolents[i]));
  }
  return result;
}

// Class to facilitate numerical integration
// 
// @param shapeH, scaleH, shapeP, scaleP double The parameters of the Weibull
//   distributions.
// @param U double The upper boundary.
//
// Initialization sets shapeH, scaleH, shapeP, scaleP, and U for each integral.
// operator() defines the integral
class WeibPDF: public Numer::Func {
private: 
    double shapeH;
    double scaleH;
    double shapeP;
    double scaleP; 
    double U;
    
public: 
    WeibPDF(double shapeH_, double scaleH_, double shapeP_, double scaleP_, double U_) : 
    shapeH(shapeH_), scaleH(scaleH_), shapeP(shapeP_), scaleP(scaleP_), U(U_){}
    double operator()(const double& x) const {
        NumericVector local_vec(1, x);
        // f(x; shape_H, scale_H) F(U - x; shape_P, scale_P)
        return dweibull(local_vec, shapeH, scaleH, false)[0] * 
            pweibull(U - local_vec, shapeP, scaleP, false, false)[0];
    } 
};

// Provided the lower and upper boundaries, integrate 
//   f(x; k, lambda) F(U - x; k, lambda)
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param L double The lower limit of integration.
// @param U double The upper limit of integration.
//
// @returns double The integral.
// @keywords internal
double compute_integral(List theta, double L, double U, double scale_P) {
    
    // Initialize instance of integral definition for the current upper bound
    WeibPDF f(theta["shape_H"], theta["scale_H"], theta["shape_P"], 
              scale_P, U);
    
    // error estimate and error code returned from the integration routine.
    // only the error flag is of interest here.
    double err_est;
    int err_code;
    
    // integrate
    double res = integrate(f, L, U, err_est, err_code);
    
    // if error code is not 0, there was a problem. abort.
    if (err_code > 0) stop("Unable to perform integration");
    
    return res;
}

// Compute integral at each age at first screening
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Ages at first screening
// @param t0 double The initial time.
//
// @returns An Rcpp::NumericVector of the integrals at each ...
// @keywords internal
NumericVector compute_cp_log(List theta, 
                             NumericVector AFS, 
                             double t0, 
                             int irateP) {
  
  double L = 0.0;        // lower bound
  NumericVector U = AFS - t0; // upper bound
    
  // F(S_0 - t_0; shape, scale)
  NumericVector prob_onset_after = pweibull(U, 
                                            theta["shape_H"], 
                                            theta["scale_H"], 
                                            false, false);
  NumericVector scale_P = theta["scale_P"];

  NumericVector integral = no_init(U.size());
  std::transform(U.begin(),
                 U.end(),
                 integral.begin(), [=] (double x) {
                   return compute_integral(theta, L, x, scale_P[irateP]); });

  double psi = theta["psi"];
  
  NumericVector cp = psi + (1.0 - psi) * (prob_onset_after + integral);
  return log(cp);
}

// Sum of the integrals
//
// @param data_object List The data for a single censor type and group.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 double The initial time.
//
// @returns double The sum of the integrals over all AFS.
//
// @keywords internal
double dloglik_cp_obj(List data_object,
                      List theta, 
                      double t0) {
  NumericMatrix age_entry = data_object["age_entry"];
  NumericVector AFS = age_entry( _ , 0);
  NumericVector nAFS = age_entry( _ , 1);
  return -sum(compute_cp_log(theta, AFS, t0, data_object["irateP"]) * nAFS);
}

// returns for each unique endpoint type/beta/rate_P combination
NumericVector dloglik_cp_List(List data_objects, 
                              List theta, 
                              double t0) {

  NumericVector result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dloglik_cp_obj(data_objects[i], theta, t0);
  }
  return result;
}

// Sum of the log-likelihood terms involving psi, the probability of
//   indolence
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 double The initial time.
//
// @returns double 
// Unk +
// sum i = 1, n
//   Grp1 ln[psi]
//   Grp2 {I(ind_i == 1) ln[psi] + 
//         I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//   Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
double dloglik_psi(List data_objects, 
                   List indolents,
                   List age_at_tau_hp_hats,
                   List theta, 
                   double t0) {
    
  // this is a vector of size data_objects
  double dlog_cp = sum(dloglik_cp_List(data_objects, theta, t0));
    
  // this is a vector of size data_objects
  double dlog_PI = sum(dloglik_PI_sum(data_objects, indolents, 
                                      age_at_tau_hp_hats, theta));
    
  return dlog_cp + dlog_PI;
}

// The log-likelihood terms involving tau for a single censor type and
//   parameter groups
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type
//   under consideration.
// @param indolent IntegerVector =1 if indolent.
// @param t0 double Initial time.
//
// @returns NumericVector(data_object["n"]) Each element provides
//   Grp1 ln[ F_h(T_i - t0) ] +
//   Grp2 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]} +
//   Grp3 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]} +
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
NumericVector dloglik_tau_obj(List data_object,
                              List theta, 
                              NumericVector age_at_tau_hp_hat, 
                              IntegerVector indolent, 
                              double t0) {
    
  // Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
  // NumericVector of length data_object["n"]
  NumericVector dlog_H = dloglik_sojourn_H_obj(data_object, theta, 
                                               age_at_tau_hp_hat, t0);

  // Grp1 {0}
  // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
  // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
  // NumericVector of length data_object["n"]
  NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta, 
                                               age_at_tau_hp_hat, indolent);
    
  // n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
  // NumericVector of length data_object["n"]
  NumericVector dlog_S = dloglik_screens_obj(data_object, theta, 
                                             age_at_tau_hp_hat);
    
  return dlog_H + dlog_P + dlog_S;
}

// The log-likelihood terms involving tau for all censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List 
//   Each element is a NumericVector providing
//   Grp1 ln[ F_h(T_i - t0) ] +
//   Grp2 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]} +
//   Grp3 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]} +
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//  for a unique endpoint type/beta/rate_P combination
//
// @keywords internal
List dloglik_tau_List(List data_objects, 
                      List indolents, 
                      List age_at_tau_hp_hats, 
                      List theta, 
                      double t0) {
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
      result[i] = dloglik_tau_obj(data_objects[i], theta, age_at_tau_hp_hats[i], 
                                  indolents[i], t0);
  }
  return result;
}

// Sum of the log-likelihood involving the rate of tau_hp distribution
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 double The initial time.
//
// @returns double
// Unk +
// sum_i = 1, n
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
double dloglik_rate_H(List data_objects, 
                      List age_at_tau_hp_hats, 
                      List theta, 
                      double t0) {
  
  // vector of size data_objects
  double dlog_cp = sum(dloglik_cp_List(data_objects, theta, t0));
    
  double dlog_H = dloglik_sojourn_H_sum(data_objects, age_at_tau_hp_hats, 
                                        theta, t0);
    
  return dlog_cp + dlog_H;
}


// The log-likelihood terms involving the rate of tau_p distribution
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 double The initial time.
//
// @returns NumericVector The result for each rate_P group
//   Unk +
//   Grp1 {0}
//   Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
//   Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
NumericVector dloglik_rate_P(List data_objects, 
                             List indolents, 
                             List age_at_tau_hp_hats,
                             List theta, 
                             double t0) {
  
  IntegerVector irateP = theta["irateP"];
  NumericVector rate_P = theta["rate_P"];
  
  NumericVector dlog_cp(rate_P.size(), 0.0);
  NumericVector dlog_P(rate_P.size(), 0.0);
  for (int i = 0; i < data_objects.size(); ++i) {
    dlog_cp[irateP[i]] += dloglik_cp_obj(data_objects[i], theta, t0);
    dlog_P[irateP[i]] += sum(dloglik_sojourn_P_obj(data_objects[i], theta, 
                                                   age_at_tau_hp_hats[i],
                                                   indolents[i]));
  }
  
  return dlog_cp + dlog_P;
  
}

//
// M-H rate_H ////////

// Sample a new rate value for the tau_hp distribution
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_H NumericVector Small value shifts.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
//
// @keywords internal
double rprop_rate_H(List theta, double epsilon_rate_H) {
    double rate_H = theta["rate_H"];
    rate_H += runif(1, - epsilon_rate_H, epsilon_rate_H)[0];
    return std::fabs(rate_H);
}

// Metropolis-Hastings step for rate of tau_HP distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_H double Small shift values. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
//
// @keywords internal
List MH_rate_H(List data_objects,
               List indolents, 
               List prior, 
               List age_at_tau_hp_hats,
               List theta_cur, 
               double epsilon_rate_H, 
               double t0) {
    
  // current value
  double rate_H_cur = theta_cur["rate_H"];
  double rate_H_new = rprop_rate_H(theta_cur, epsilon_rate_H);
  NumericVector rates = NumericVector::create(rate_H_new, rate_H_cur);

  List theta_new = clone(theta_cur);
  theta_new = add_rate_H(theta_new, rate_H_new);  // for dloglik_rate_H()
    
  // Unk +
  // sum_i = 1, n
  //   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
  // calculated with current rate_H value
  double dlog_lik_cur = dloglik_rate_H(data_objects, age_at_tau_hp_hats, 
                                       theta_cur, t0);
  // Unk +
  // sum_i = 1, n
  //   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
  // calculated with new rate_H value
  double dlog_lik_new = dloglik_rate_H(data_objects, age_at_tau_hp_hats, 
                                       theta_new, t0);
    
  double prior_rate = prior["rate_H"];
    
  NumericVector dlog_prior = dgamma(rates, prior["shape_H"], 1.0 / prior_rate, true);

  double MH_logratio = (dlog_lik_new + dlog_prior[0]) - 
                       (dlog_lik_cur + dlog_prior[1]);

  double accept_prob = std::exp(MH_logratio);
  
  if (runif(1)[0] < accept_prob) {
    return List::create(Named("theta") = theta_new, 
                        Named("accept") = true,
                        Named("probability") = accept_prob);
  } else {
    return List::create(Named("theta") = theta_cur, 
                        Named("accept") = false,
                        Named("probability") = accept_prob);
  }
}

//
// M-H rate_P ////////

// Sample a new rate value for the tau_p distribution
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_P NumericVector A small value shift for each rate_p parameter.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
//
// @keywords internal
NumericVector rprop_rate_P(List theta, NumericVector epsilon_rate_P) {
  NumericVector rate_P = theta["rate_P"];
  NumericVector rate_P_new = clone(rate_P);
  for (int i = 0; i < epsilon_rate_P.size(); ++i) {
    rate_P_new[i] += runif(1, -epsilon_rate_P[i], epsilon_rate_P[i])[0];
  }
  return abs(rate_P_new);
}

// Metropolis-Hastings step for rate of tau_P distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements.  
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_P NumericVector A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
//
// @keywords internal
List MH_rate_P(List data_objects,
               List indolents, 
               List prior, 
               List age_at_tau_hp_hats,
               List theta_cur, 
               NumericVector epsilon_rate_P, 
               double t0) {

  NumericVector rate_P_cur = theta_cur["rate_P"];
  NumericVector rate_P_new = rprop_rate_P(theta_cur, epsilon_rate_P); 

  List theta_new = clone(theta_cur);
  theta_new = add_rate_P(theta_new, rate_P_new); 
  
  // M-H acceptance ratio
  
  // Unk +
  // Grp1 {0}
  // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
  // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
  // evaluated at current rate_P
  // NumericVector of lengthe rate_P
  NumericVector dlog_lik_cur = dloglik_rate_P(data_objects, indolents, 
                                              age_at_tau_hp_hats, 
                                              theta_cur, t0);
  // Unk +
  // Grp1 {0}
  // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
  // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
  // evaluated at new rate_P
  // NumericVector of lengthe rate_P
  NumericVector dlog_lik_new = dloglik_rate_P(data_objects, indolents, 
                                              age_at_tau_hp_hats,
                                              theta_new, t0);
  
  NumericVector prior_rate = prior["rate_P"];
  NumericVector prior_shape = prior["shape_P"];
  
  NumericVector dlog_prior_new(rate_P_new.size(), 0.0);
  NumericVector dlog_prior_cur(rate_P_cur.size(), 0.0);
  NumericVector tmp(1);
  for (int i = 0; i < rate_P_cur.size(); ++i) {
    tmp[0] = rate_P_new[i];
    dlog_prior_new[i] = dgamma(tmp, prior_shape[i], 
                               1.0 / prior_rate[i], true)[0];
    tmp[0] = rate_P_cur[i];
    dlog_prior_cur[i] = dgamma(tmp, prior_shape[i], 
                               1.0 / prior_rate[i], true)[0];
  }

  NumericVector MH_logratio = (dlog_lik_new + dlog_prior_new) - 
                              (dlog_lik_cur + dlog_prior_cur);

  NumericVector accept_prob = exp(MH_logratio);
  LogicalVector test = runif(MH_logratio.size()) < accept_prob;
  
  NumericVector accepted_rate_P = rate_P_cur;
  accepted_rate_P[test] = rate_P_new[test];
  
  List theta_accepted = clone(theta_cur);
  theta_accepted = add_rate_P(theta_accepted, accepted_rate_P);

  return List::create(Named("theta") = theta_accepted, 
                      Named("accept") = test,
                      Named("probability") = accept_prob);
}


//
// M-H age_at_tau_hp_hat ////////

// Compute the probability of tau falling within the observed intervals for a 
//   single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @return List Elements include `values`, a NumericVector of all probabilities,
//   `starts`, an IntegerVector containing the first element of `values
//   pertaining to case i; `ends`, an IntegerVector containing the last element
//   pertaining to individual i; and `lengths` the number of element in `values`
//   pertaining to individual i.
//
// @keywords internal
// [[Rcpp::export]]
List compute_prob_tau_obj(List data_object, List theta, double t0) {
  
  NumericVector endpoint_time = data_object["endpoint_time"];

  // extract vectorized observed interval data
  // note that these interval boundaries are actually "age at time of"
  // the interval is obtained by subtracting the lower age if no clinical
  // diagnosis or by subtracting from the age at clinical diagnosis
  List endpoints = data_object["endpoints"];

  NumericVector ep = endpoints["values"];
  IntegerVector starts = endpoints["starts"];
  IntegerVector ends = endpoints["ends"];
  IntegerVector lengths = endpoints["lengths"];

  
  int endpoint_type = data_object["endpoint_type"];

  // estimating the probability of each interval
  // The number of intervals is 1 less than the number of boundary points
  IntegerVector pt_lengths = lengths - 1;
  // The final index for each individual is 1 fewer than the lengths (0 indexing)
  IntegerVector pt_ends = cumsum(pt_lengths);
  pt_ends = pt_ends - 1;
  // The starts are obtained using the final index and the lengths
  IntegerVector pt_starts = pt_ends - pt_lengths + 1;
  
  // compute probability of interval based on the Weibull waiting time in H 
  // (screen, censored) or in P (clinical)
  NumericVector prob_interval;
  
  // paired points at which the distribution is evaluated 
  NumericVector vec1(sum(pt_lengths));
  NumericVector vec2(sum(pt_lengths));
  int j;
  int cnt = 0;
  
  for (int i = 0; i < starts.length(); ++i) {
    j = starts[i];
    while (j < ends[i]) {
      vec1[cnt] = ep[j];
      vec2[cnt] = ep[j+1];
      cnt++;
      j++;
    }
  }

  if (endpoint_type != 3) {
    // Preclinical and censored cases
    // Pr(interval_k_lower < tau_hp <= interval_k_upper)
    prob_interval = pweibull_ab(vec1 - t0, vec2 - t0, 
                                theta["shape_H"], theta["scale_H"]);
    
  } else {
    // Clinical cases
    int irateP = data_object["irateP"];
    NumericVector scale_P = theta["scale_P"];
    
    // need to repeat censor time to match interval data structure
    NumericVector ct = repVec(endpoint_time, pt_lengths);
    // Pr(interval_l_lower < tau_p <= interval_l_upper)
    prob_interval = pweibull_ab(ct - vec2,
                                ct - vec1,
                                theta["shape_P"], scale_P[irateP]);
  }

  NumericVector prob_tau = prob_interval;
  NumericVector norm;
  
  // we need to repeat (1-beta)^{K-1}:0
  List ages_screen = data_object["ages_screen"];
  if (ages_screen.size() == 0) {
    // if not screening history, normalize the probabilities
    // without modification by screen sensitivities
    for (int i = 0; i < pt_starts.size(); ++i) {
      norm = prob_tau[seq(pt_starts[i], pt_ends[i])];
      norm = norm / sum(norm);
      prob_tau[seq(pt_starts[i], pt_ends[i])] = norm;
    }
    
    return List::create(Named("values") = prob_tau,
                        Named("starts") = pt_starts,
                        Named("ends") = pt_ends,
                        Named("lengths") = pt_lengths);
    
  }
  
  IntegerVector screen_types = ages_screen["types"];
  IntegerVector age_starts = ages_screen["starts"];
  IntegerVector age_ends = ages_screen["ends"];
  IntegerVector age_lengths = ages_screen["lengths"];
  NumericVector betas = theta["beta"];

  double res;
  int tau_index = 0;
  for (int i = 0; i < pt_starts.size(); ++i) {
      
    if (pt_starts[i] == pt_ends[i]) {
      // if only 1 interval set to 1
      prob_tau[pt_starts[i]] = 1.0;
      continue; 
    }
    
    // start from the final interval and final screen age
    tau_index = pt_ends[i];
    j = age_ends[i];
    
    if (pt_lengths[i] == age_lengths[i]) {
      // means that the endpoint time was a screen
      // ep = t0, S1, S2, S3, TE -- TE == S3 for screen detected
      // ep = t0, S1, S2, S3
      // intervals are [t0, S1], [S1, S2], [S2, S3]
      res = betas[screen_types[j]];
      prob_tau[tau_index] *= res;
      tau_index--;
      j--;
      while (j >= age_starts[i]) {
        res *= (1.0 - betas[screen_types[j]]);
        prob_tau[tau_index] *= res;
        tau_index--;
        j--;
      }
    } else {
      res = 1.0;
      while (j >= age_starts[i]) {
        prob_tau[tau_index] *= res;
        res *= (1.0 - betas[screen_types[j]]);
        tau_index--;
        j--;
      }
      prob_tau[pt_starts[i]] *= res;
    }
    
    norm = prob_tau[seq(pt_starts[i], pt_ends[i])];
    norm = norm / sum(norm);
    
    prob_tau[seq(pt_starts[i], pt_ends[i])] = norm;
  }
  
  return List::create(Named("values") = prob_tau,
                      Named("starts") = pt_starts,
                      Named("ends") = pt_ends,
                      Named("lengths") = pt_lengths);
}

// Compute probability of tau for all censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is itself a List
//   containing the probabilities of tau details for the specific
//   censoring type.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List compute_prob_tau_List(List data_objects, List theta, double t0) {
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = compute_prob_tau_obj(data_objects[i], theta, t0);
  }
  return result;
}

double compute_optimal_lambda(double shape_P, double rate_P) {
  
  if (shape_P < 1.0001) {
    return rate_P;
  }
  
  double d = shape_P / (shape_P - 1.0);
  double tmp = 1.0 / (1.0 - shape_P);
  double b1 = - pow(shape_P, shape_P * tmp) * pow(rate_P, tmp);
  double b2 = pow(shape_P, tmp) * pow(rate_P, tmp);

  return pow((b1 + b2) * d, -1.0 / d);
    
}

// Sample current estimated probability of tau for new tau_hp
//  for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns NumericVector The new estimated tau for each case in the 
//   censoring type under consideration
//
// @keywords internal
// [[Rcpp::export]]
NumericVector rprop_age_at_tau_hp_hat_obj(List data_object, 
                                          List prob_tau, 
                                          List theta, 
                                          double t0) {
    
  NumericVector endpoint_time = data_object["endpoint_time"];
    
  // extract vectorized interval information from data list
  // note that these are "age at time of" boundaries
  List endpoints = data_object["endpoints"];
  NumericVector ep = endpoints["values"];
  IntegerVector ep_starts = endpoints["starts"];
  IntegerVector ep_ends = endpoints["ends"];
  IntegerVector ep_lengths = endpoints["lengths"];
    
  int endpoint_type = data_object["endpoint_type"];
    
  // extract vectorized probabilities of tau
  NumericVector pt = prob_tau["values"];
  IntegerVector pt_starts = prob_tau["starts"];
  IntegerVector pt_ends = prob_tau["ends"];
  IntegerVector pt_lengths = prob_tau["lengths"];
    
  // sample each case's interval
    
  IntegerVector k_new = no_init(endpoint_time.length());
  for (R_xlen_t i = 0; i < pt_starts.length(); ++i) {
    int K = pt_lengths[i];
    NumericVector vec = pt[seq(pt_starts[i], pt_ends[i])];
    k_new[i] = sample(K, 1, false, vec)[0] - 1;
  }
    
  NumericVector age_at_tau_hp_hat_new;

  NumericVector vec1 = ep[ep_starts + k_new];
  NumericVector vec2 = ep[ep_starts + k_new + 1];

  if (endpoint_type == 1) {
    // Preclinical cases
    // sample tau_hp from (interval_lower, interval_upper]
    NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                 theta["shape_H"], 
                                                 theta["scale_H"]);
    age_at_tau_hp_hat_new = t0 + sojourn_H_new;
  } else if (endpoint_type == 2) {
    // Censored w/ screening data cases
    // sample tau_hp from (interval_lower, interval_upper]
    NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                 theta["shape_H"], 
                                                 theta["scale_H"]);
    age_at_tau_hp_hat_new = t0 + sojourn_H_new;
        
    // if the estimated age at tau_hp > endpoint_time set tau = Inf
    LogicalVector tst = age_at_tau_hp_hat_new >= endpoint_time;
    age_at_tau_hp_hat_new[tst] = R_PosInf;
  } else if (endpoint_type == 3) {
    // Clinical cases
    int irateP = data_object["irateP"];
    NumericVector scale_P = theta["scale_P"];
    
    // sample tau_p from (interval_lower, interval_upper]
    NumericVector sojourn_P_new = rweibull_trunc(endpoint_time - vec2, 
                                                 endpoint_time - vec1,
                                                 theta["shape_P"], 
                                                 scale_P[irateP]);
    // age at tau_hp is age at clinical diagnosis - tau_p
    age_at_tau_hp_hat_new = endpoint_time - sojourn_P_new;
  } else if (endpoint_type == 4) {
    // Censored cases that have no screens
    // for these individuals, endpoint = (t0, T, Inf)
    // for individuals with k_new = 1, the sampling interval is in the
    // |T, inf) interval. These are unchanged.
    // for individuals with k_new = 0, the sampling interval is in the
    // |t0, T) interval. These will now sample from an exponential
    
    int n = data_object["n"];
    age_at_tau_hp_hat_new = rep(R_PosInf, n);
    
    LogicalVector finite_intv = k_new == 0;
    if (is_true(any(finite_intv))) {
      int irateP = data_object["irateP"];
      NumericVector rate_P = theta["rate_P"];
      double lambda = compute_optimal_lambda(theta["shape_P"], rate_P[irateP]);
      int n_finite = sum(as<IntegerVector>(finite_intv));
      vec2 = rexp(n_finite, lambda);
      vec1 = endpoint_time[finite_intv];
      NumericVector vec3 = vec1 - vec2;
      age_at_tau_hp_hat_new[finite_intv] = vec3;
    }

  }
  return age_at_tau_hp_hat_new;
}


// Sample current estimated probability of tau for new tau_hp for all censor
//  censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List rprop_age_at_tau_hp_hat_List(List data_objects, 
                                  List prob_tau, 
                                  List theta, 
                                  double t0) {
  
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = rprop_age_at_tau_hp_hat_obj(data_objects[i], prob_tau[i], 
                                            theta, t0);
  }
  return result;
}

// Contribution to the log-likelihood from probability of tau
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical for all participants of the censor type under 
//   consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector 
//
// @keywords internal
NumericVector dlog_prop_age_at_tau_hp_hat_obj(List data_object,
                                              List prob_tau,
                                              List theta,
                                              NumericVector age_at_tau_hp_hat, 
                                              double t0) {
    
  NumericVector endpoint_time = data_object["endpoint_time"];

  // extract vectorized interval boundary information
  // note that this an "age at time of" vector
  List endpoints = data_object["endpoints"];
  NumericVector ep = endpoints["values"];
  IntegerVector ep_starts = endpoints["starts"];
  IntegerVector ep_ends = endpoints["ends"];
  IntegerVector ep_lengths = endpoints["lengths"];
    
  int endpoint_type = data_object["endpoint_type"];

  // extract vectorized probability of tau information
  NumericVector pt = prob_tau["values"];
  IntegerVector pt_starts = prob_tau["starts"];
  IntegerVector pt_ends = prob_tau["ends"];
  IntegerVector pt_lengths = prob_tau["lengths"];
  
  NumericVector dlog_tau;
  NumericVector dlog_k;
  if (endpoint_type != 4) {
    
    // replicate age_at_tau_hp_hat_i to correlate with vectorized boundaries
    NumericVector tmp_age_at_tau = repVec(age_at_tau_hp_hat, ep_lengths);
    
    LogicalVector cmp = ep < tmp_age_at_tau;
    IntegerVector cmp_csum = cumsum(as<IntegerVector>(cmp));
    LogicalVector cmp_starts = cmp[ep_starts];
    IntegerVector icmp_starts = as<IntegerVector>(cmp_starts);
    
    // k_new is the last screen before the current estimated transition time
    IntegerVector k_new = cmp_csum[ep_ends] - cmp_csum[ep_starts] + icmp_starts - 1;
  
    LogicalVector test = ((ep_starts + k_new + 1) >= ep.size()) |
      (k_new < 0);
    if (is_true(any(test))) {
      stop("k_new is going out of bounds; contact developer");
    }
  
    test = ((pt_starts + k_new) >= pt.size()) |
      ((pt_starts + k_new) < 0);
    if (is_true(any(test))) {
      stop("k_new is going out of boundspf pt; contact developer");
    } 
  
    // extract the corresponding probability of tau
    NumericVector tmp_vec = pt[pt_starts + k_new];
    dlog_k = log(tmp_vec);
  
    NumericVector vec1 = ep[ep_starts + k_new];
    NumericVector vec2 = ep[ep_starts + k_new + 1];
    
    if (endpoint_type == 1) {
      // Preclinical cases
      // P(tau_hp_i == tau_hp_hat_i | interval_lower < tau_hp <= interval_upper)
      dlog_tau = dweibull_trunc(age_at_tau_hp_hat - t0, 
                                vec1 - t0, 
                                vec2 - t0,
                                theta["shape_H"], theta["scale_H"], true);
    } else if (endpoint_type == 2) {
      // Censored cases w/ screens
      // P(tau_hp_i == tau_hp_hat_i | interval_lower < tau_hp <= interval_upper)
      dlog_tau = dweibull_trunc(age_at_tau_hp_hat - t0, 
                                vec1 - t0, 
                                vec2 - t0,
                                theta["shape_H"], theta["scale_H"], true);
      // if estimated age at tau_hp > censor time, set to 0
      LogicalVector tst = age_at_tau_hp_hat >= endpoint_time;
      dlog_tau[tst] = 0.0;
    } else if (endpoint_type == 3) {
      // Clinical cases
      int irateP = data_object["irateP"];
      NumericVector scale_P = theta["scale_P"];

      // P(tau_p_i == tau_p_hat_i | interval_lower < tau_p <= interval_upper)
      dlog_tau = dweibull_trunc(endpoint_time - age_at_tau_hp_hat, 
                                endpoint_time - vec2, 
                                endpoint_time - vec1,
                                theta["shape_P"], scale_P[irateP], true);
    }
  } else if (endpoint_type == 4) {
    // Censored cases w/out screens
    int n = data_object["n"];
    
    dlog_tau = rep(0.0, n);
    NumericVector tmp_vec = pt[pt_ends];
    
    LogicalVector tst = age_at_tau_hp_hat < endpoint_time;
    
    if (is_true(any(tst))) {
      int irateP = data_object["irateP"];
      NumericVector rate_P = theta["rate_P"];
      double lambda = compute_optimal_lambda(theta["shape_P"], rate_P[irateP]);
      NumericVector vec1 = endpoint_time[tst] - age_at_tau_hp_hat[tst];
      NumericVector vec2 = dexp(vec1, lambda, true);

      dlog_tau[tst] = vec2;
      NumericVector start_pt = pt[pt_starts];
      tmp_vec[tst] = start_pt[tst];
      
    }    
    dlog_k = log(tmp_vec);
    
  } else {
    stop("unrecognized endpoint type");
  }
  
  return dlog_k + dlog_tau;
}

// The log-likelihood terms involving probability of tau for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped
//   according to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// @keywords internal
List dlog_prop_age_at_tau_hp_hat_List(List data_objects, 
                                      List prob_taus, 
                                      List age_at_tau_hp_hats, 
                                      List theta, double t0) {
  
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dlog_prop_age_at_tau_hp_hat_obj(data_objects[i], prob_taus[i], 
                                                theta, age_at_tau_hp_hats[i], 
                                                t0);
  }
  return result;
}

//
// indolent ////////

// Probability of being indolent for each case of a single censor type given 
//   current parameter and transition time estimates
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type
//   under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector 
//
// @keywords internal
// [[Rcpp::export]]
NumericVector compute_prob_indolent_obj(List data_object, List theta,
                                        NumericVector age_at_tau_hp_hat) {
    
  IntegerVector ind(age_at_tau_hp_hat.size(), 0);
    
  // Grp1 ln[psi]
  // Grp2 {I(ind_i == 1) ln[psi] + 
  //       I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
  // Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
  // evaluated when all cases are set as non-indolent
  // NumericVector of length data_object["n"]
  NumericVector L_0 = dloglik_PI_obj(data_object, theta, age_at_tau_hp_hat, ind);
  L_0 = exp(L_0);
    
  std::fill(ind.begin(), ind.end(), 1);
  // Grp1 ln[psi]
  // Grp2 {I(ind_i == 1) ln[psi] + 
  //       I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
  // Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
  // evaluated when all cases are set as indolent
  // NumericVector of length data_object["n"]
  NumericVector L_1 = dloglik_PI_obj(data_object, theta, age_at_tau_hp_hat, ind);
  L_1 = exp(L_1);
    
  return L_1 / (L_0 + L_1);
}

// Probability of being indolent for all cases given current parameter and
//  transition time estimates
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   according to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List compute_prob_indolent_List(List data_objects, List age_at_tau_hp_hats, List theta) {
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
      result[i] = compute_prob_indolent_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
  }
  return result;
}

// Generate vector of indolence for all cases of a single censor type based on
//   current estimated probability of indolence.
//
// Note: indolence is always 0 (false) for clinical cases
//
// @param data_object List The data for a single censoring type.
// @param prob_indolent NumericVector the probability of each case being
//   indolent.
//
// @returns IntegerVector Indicator of indolence.
//
// @keywords internal
// [[Rcpp::export]]
IntegerVector rprop_indolent_obj(List data_object, NumericVector prob_indolent) {
    
  int endpoint_type = data_object["endpoint_type"];
    
  int n = prob_indolent.size();
  IntegerVector indolent(n);
    
  if (endpoint_type == 3) {
    // Clinical cases
    // indolence is always false for those with a clinical diagnosis
    indolent.fill(0);
  } else {
    // Censored and preclinical cases
    // indolence is a random sample from binomial using the current
    // estimated probability of indolence
    std::transform(prob_indolent.begin(), 
                   prob_indolent.end(), 
                   indolent.begin(), [=](double p){ return rbinom(1, 1, p)[0]; }); 
  }
    
  return indolent;
}

// Generate vector of indolence for all cases.
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_indolents List The current indolence probability broken down
//   according to the censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
//
// @returns List There are 3 elements. Each element is an IntegerVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List rprop_indolent_List(List data_objects, List prob_indolents) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = rprop_indolent_obj(data_objects[i], prob_indolents[i]);
    }
    return result;
}


// Case contributions to the log-likelihood for the probability of indolence
//   for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param prob_indolent NumericVector the probability of each case being
//   indolent.
// @param indolent IntegerVector the current indolence indicator
//
// @returns NumericVector 
//
// {Grp1 + Grp2} {I(ind_i == 1) ln[ Pr(Ind == 1) ] + 
//                I(ind_i == 0) ln[ 1 - Pr(Ind == 1) ]} + Grp3 {0}
//
// @keywords internal
// [[Rcpp::export]]
NumericVector dlog_prop_indolent_obj(List data_object, 
                                     NumericVector prob_indolent,
                                     IntegerVector indolent) {
    
  int endpoint_type = data_object["endpoint_type"];
    
  NumericVector result = log(1.0 - prob_indolent);
  if (endpoint_type == 3) {
    // Clinical cases
    // dlog-likelihood is 0 for all clinical cases
    result.fill(0.0);
  } else {
    // Censored and preclinical cases
    // dlog-likelihood is log(1-p) for 0 cases; log(p) for 1 cases
    LogicalVector ind1 = indolent == 1;
    NumericVector prob1 = prob_indolent[ind1];
    prob1 = log(prob1);
    result[ind1] = prob1;
  }
    
  return result;
}

// Log-likelihood for the probability of indolence for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current indolence indicator broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_indolents List The current probability of indolence broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
//
// @returns double
// sum i = 1, n
//   {Grp1 + Grp2} {I(ind_i == 1) ln[ Pr(Ind == 1) ] + 
//                  I(ind_i == 0) ln[ 1 - Pr(Ind == 1) ]} + Grp3 {0}
//
// @keywords internal
double dlog_prop_indolent_sum(List data_objects, 
                              List indolents,
                              List prob_indolents) {
  double result = 0.0;
  for (int i = 0; i < data_objects.size(); ++i) {
    result += sum(dlog_prop_indolent_obj(data_objects[i], prob_indolents[i], 
                                         indolents[i]));
  }
  return result;
}

// Not used
// [[Rcpp::export]]
NumericVector dlog_prop_latent_obj(List data_object, 
                                   List prob_tau,
                                   List theta,
                                   NumericVector age_at_tau_hp_hat,
                                   NumericVector prob_indolent,
                                   IntegerVector indolent, 
                                   double t0) {
  
  // vector of length data_object["n"]
  NumericVector dlog_prop_tau = dlog_prop_age_at_tau_hp_hat_obj(data_object,
                                                                prob_tau,
                                                                theta,
                                                                age_at_tau_hp_hat, 
                                                                t0);
  
  // vector of length data_object["n"]
  NumericVector dlog_prop_I = dlog_prop_indolent_obj(data_object, 
                                                     prob_indolent,
                                                     indolent);
  
  return dlog_prop_tau + dlog_prop_I;
  
}

// Not used
List dlog_prop_latent_List(List data_objects, 
                           List indolents, 
                           List age_at_tau_hp_hats, 
                           List prob_indolents,
                           List prob_taus,
                           List theta,
                           double t0) {
  
  List result(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    result[i] = dlog_prop_latent_obj(data_objects[i], 
                                     prob_taus[i],
                                     theta,
                                     age_at_tau_hp_hats[i],
                                     prob_indolents[i],
                                     indolents[i], 
                                     t0);
  }
  
  return result;  
}

//
// psi ////////

// Random walk in psi space
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_psi A double. The maximum step size.
//
// @returns A double. The new psi value
//
// @keywords internal
double rprop_psi(List theta, double epsilon_psi) {
    
  double psi = theta["psi"];
    
  // uniform random walk
  double psi_prop = runif(1, psi - epsilon_psi, psi + epsilon_psi)[0];
  while(psi_prop < 0.0 || psi_prop > 1.0) {
    // reflection on lower bound 0 and upper bound 1 
    if (psi_prop < 0.0) {
      psi_prop = 0.0 + (0.0 - psi_prop);
    } else if (psi_prop > 1.0) {
      psi_prop = 1.0 - (psi_prop - 1.0);
    }
  }
    
  return psi_prop;
}

//
// M-H tau ////////

// Metropolis-Hastings step for tau distributions for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector Indicator of indolence.
// @param t0 A scalar double The initial time.
//
// @returns List Element `$age_at_tau_hp_hat` is the possibly updated List of
//   tau values and element `$accepted` is a LogicalVector, where true indicates
//   that the corresponding element of `$age_at_tau_hp_hat` has been updated.
//
// @keywords internal
List MH_tau_obj(List data_object, List theta, NumericVector age_at_tau_hp_hat, 
                IntegerVector indolent, double t0) {
 
  // propose new latent variables age_at_tau_hp_hat/tau_CP
  List prob_tau = compute_prob_tau_obj(data_object, theta, t0);

  NumericVector age_at_tau_hp_hat_new = rprop_age_at_tau_hp_hat_obj(data_object, 
                                                                    prob_tau, 
                                                                    theta, t0);

  // M-H acceptance ratio
  NumericVector dlog_prop_cur = dlog_prop_age_at_tau_hp_hat_obj(data_object,
                                                                prob_tau,
                                                                theta,
                                                                age_at_tau_hp_hat, 
                                                                t0);

  NumericVector dlog_prop_new = dlog_prop_age_at_tau_hp_hat_obj(data_object,
                                                                prob_tau,
                                                                theta,
                                                                age_at_tau_hp_hat_new,
                                                                t0);

  NumericVector dlog_lik_cur = dloglik_tau_obj(data_object,
                                               theta,
                                               age_at_tau_hp_hat, 
                                               indolent, 
                                               t0);

  NumericVector dlog_lik_new = dloglik_tau_obj(data_object,
                                               theta,
                                               age_at_tau_hp_hat_new, 
                                               indolent, 
                                               t0);

  NumericVector MH_logratio = (dlog_lik_new - dlog_prop_new) - 
        (dlog_lik_cur - dlog_prop_cur);

  LogicalVector test = runif(MH_logratio.size()) < exp(MH_logratio);
    
  NumericVector accepted_tau = age_at_tau_hp_hat;
  accepted_tau[test] = age_at_tau_hp_hat_new[test];
  
  // use NA to indicate Inf -> Inf cases
  LogicalVector na_cases = is_infinite(age_at_tau_hp_hat) & is_infinite(age_at_tau_hp_hat_new);
  test[na_cases] = NA_LOGICAL;
  
  return List::create(Named("age_at_tau_hp_hat") = accepted_tau,
                      Named("accept") = test);
    
}

// Metropolis-Hastings step for rate of tau distribution for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current indolence indicator broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical"
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List Element `$age_at_tau_hp_hats` is the updated estimates for tau 
//   broken down according to the censor type and element `$accept` is a list 
//   broken down according to the censor type indicating if the corresponding
//   elements of `$age_at_tau_hp_hats` were updated.
//
// @keywords internal
List MH_tau_List(List data_objects, 
                 List indolents, 
                 List age_at_tau_hp_hats, 
                 List theta, 
                 double t0) {
    
  List result_tau(data_objects.size());
  List result_accept(data_objects.size());
  for (int i = 0; i < data_objects.size(); ++i) {
    List res = MH_tau_obj(data_objects[i], theta, age_at_tau_hp_hats[i], 
                          indolents[i], t0);
    result_tau[i] = res["age_at_tau_hp_hat"];
    result_accept[i] = res["accept"];
  }
    
  return List::create(Named("age_at_tau_hp_hats") = result_tau,
                      Named("accept") = result_accept);  
}


//
// M-H (psi, indolent) ////////

// Sample current probability of indolence distribution for each case of a 
//   single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> preclinical transition for all participants of the censor type
//   under consideration.
//
// @returns List Element `$indolent` is the updated indolence indicator and
//   element `$dlog_prop` is the updated derivative evaluated at the updated
//   indolence indicators.
//
// @keywords internal
List rprop_dlog_indolent_obj(List data_object, List theta, NumericVector age_at_tau_hp_hat) {
    
    NumericVector prob_indolent_new = compute_prob_indolent_obj(data_object,
                                                                theta,
                                                                age_at_tau_hp_hat);
    // Generate new indolence vector
    IntegerVector indolent_new = rprop_indolent_obj(data_object, prob_indolent_new);

    double dlog_prop_new = sum(dlog_prop_indolent_obj(data_object, 
                                                      prob_indolent_new, 
                                                      indolent_new));
    
    return List::create(Named("indolent") = indolent_new,
                        Named("dlog_prop") = dlog_prop_new);
}

// Sample current probability of indolence distribution for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List Element `$indolent_new` is a List of the the updated 
//   indolence indicator for each censor type and element `$dlog_prop_indolent_new`
//   is the updated derivative for each censor type.
//
// @keywords internal
List rprop_dlog_indolent_List(List data_objects, List age_at_tau_hp_hats, List theta) {
    List result_indolent(data_objects.size());
    double dlog_prop = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        List res = rprop_dlog_indolent_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
        result_indolent[i] = res["indolent"];
        double dlog = res["dlog_prop"];
        dlog_prop += dlog;
    }
    
    return List::create(Named("indolent_new") = result_indolent,
                        Named("dlog_prop_indolent_new") = dlog_prop);  
}

// Metropolis-Hastings step for psi distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical"
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta_cur List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_psi A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 3 elements. `$indolence` the possibly updated indolence
//   indicator broken down by censor type; `$theta` the accepted new theta; and
//   `$accept`, true if theta and indolence were updated; FALSE otherwise.
//
// @keywords internal
List MH_psi_indolent(List data_objects,
                     List indolents, 
                     List prior, 
                     List age_at_tau_hp_hats,
                     List theta_cur, 
                     double epsilon_psi, 
                     double t0) {
    
  // propose psi
  double psi_new = rprop_psi(theta_cur, epsilon_psi);
  List theta_new = clone(theta_cur);
  theta_new = add_psi(theta_new, psi_new);
    
  // prior density
  double cur_psi = theta_cur["psi"];
  double a_psi = prior["a_psi"];
  double b_psi = prior["b_psi"];
  NumericVector tmp = NumericVector::create(psi_new, cur_psi);
  NumericVector dlog_prior = dbeta(tmp, a_psi, b_psi, true);

  // propose indolent and compute log proposal density
  List out = rprop_dlog_indolent_List(data_objects, age_at_tau_hp_hats, 
                                        theta_new);
  List indolent_new = out["indolent_new"];
    
  // proposal density
  double dlog_prop_indolent_new = out["dlog_prop_indolent_new"];
  
  List prob_indolent_cur = compute_prob_indolent_List(data_objects, 
                                                      age_at_tau_hp_hats, 
                                                      theta_cur);
  double dlog_prop_indolent_cur = dlog_prop_indolent_sum(data_objects,
                                                         indolents, 
                                                         prob_indolent_cur);
    
  // log likelihood
  double dlog_lik_cur = dloglik_psi(data_objects, indolents, 
                                    age_at_tau_hp_hats, 
                                    theta_cur, t0);
  double dlog_lik_new = dloglik_psi(data_objects, indolent_new, 
                                    age_at_tau_hp_hats, 
                                    theta_new, t0);
    
  // M-H acceptance ratio
  double MH_logratio = (dlog_lik_new + dlog_prior[0] - dlog_prop_indolent_new) - 
      (dlog_lik_cur + dlog_prior[1] - dlog_prop_indolent_cur);

  double accept_prob = std::exp(MH_logratio);
  if (runif(1)[0] < accept_prob) {
      return List::create(Named("indolents") = indolent_new, 
                          Named("theta") = theta_new, 
                          Named("accept") = true,
                          Named("probability") = accept_prob);
  } else {
      return List::create(Named("indolents") = indolents, 
                          Named("theta") = theta_cur, 
                          Named("accept") = false,
                          Named("probability") = accept_prob);
  }
    
}


//
// MCMC ////////

// Full Markov Chain Monte Carlo Procedure
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> preclinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param double epsilon_rate_H A positive scalar double specifying the maximum
//   step size in the rate_H space.
// @param double epsilon_rate_P A positive NumericVector specifying the maximum
//   step size in the rate_P space for each element of rate_P.
// @param double epsilon_psi A positive scalar double specifying the maximum
//   step size in the psi space.
// @param t0 double The initial time.
// @param M A positive scalar integer specifying the number of MC samples.
// @param thin A positive scalar integer specifying the number of samples to
//   be skipped between accepted samples.
// @param M_thin A positive scalar integer specifying the total number of
//   samples that will be accepted.
// @param n_obs A positive scalar integer specifying the total number of
//   participants. Note this should be equivalent to the sum of element
//   `$n` in each element of `data_objects`. If not saving latent
//   variables, this will be 1L
// @param n_screen_positive_total A positive IntegerVector indicating the
//   number of positive screens in each beta subset.
//
// @returns List 
//
// exported main function
// [[Rcpp::export]]
List MCMC_cpp_internal(List data_objects, List indolents, List prior, 
                       List age_at_tau_hp_hats, List theta,
                       double epsilon_rate_H, NumericVector epsilon_rate_P, double epsilon_psi,
                       double t0, int M, int thin, int M_thin, int n_obs,
                       IntegerVector n_screen_positive_total,
                       List adaptive,
                       int verbose, int save_latent) {

  NumericVector beta = theta["beta"];
  int n_rateP = epsilon_rate_P.size();
  int n_beta = beta.size();
  
  int latent_var_size = M_thin;
  if (save_latent == 0) latent_var_size = 1;
  
  NumericVector kRATE_H(M_thin), kPSI(M_thin);
  LogicalVector kACCEPT_PSI(M_thin), kACCEPT_RATE_H(M_thin);
  NumericMatrix kage_at_tau_hp_hat(latent_var_size, n_obs);
  NumericMatrix kBETA(M_thin, n_beta), kRATE_P(M_thin, n_rateP);
  IntegerMatrix kINDOLENT(latent_var_size, n_obs);
  LogicalMatrix kACCEPT_RATE_P(M_thin, n_rateP);
  LogicalMatrix kACCEPT_LATENT(latent_var_size, n_obs);
  int ikept = -1;
  
  // Used to re-sort kept data
  NumericVector tmp;
  int total_cases = 0;
  IntegerVector cases(age_at_tau_hp_hats.size(), 0);
  for(R_xlen_t k = 0; k < age_at_tau_hp_hats.length(); ++k) {
    tmp = age_at_tau_hp_hats[k];
    cases[k] = tmp.length();
    total_cases += cases[k];
  }
  
  NumericVector age_at_tau_hp_hat_thin(total_cases);
  IntegerVector indolent_thin(total_cases);
  LogicalVector accept_latent_thin(total_cases);
  
  double delta = adaptive["delta"];
  int warmup = adaptive["warmup"];
  double kappa = adaptive["kappa"];
  double gamma = adaptive["gamma"];
  int m0 = adaptive["m0"];
  
  double H_psi = 0.0;
  double H_rateH = 0.0;
  NumericVector H_rateP(epsilon_rate_P.size(), 0.0);

  int m_mod = M;
  if (verbose == 1) {
    m_mod = floor(M * 0.02);
    if (m_mod > 0) {
      Rcout << "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << std::endl;
      Rcout << "[----|----|----|----|----|----|----|----|----|----|" << std::endl;
    }
  }
  
  NumericVector adaptive_epsilon_psi;
  NumericVector adaptive_epsilon_rateH;
  NumericMatrix adaptive_epsilon_rateP;
  if (warmup > 0) {
    adaptive_epsilon_psi = NumericVector(warmup);
    adaptive_epsilon_rateP = NumericMatrix(warmup, epsilon_rate_P.size());
    adaptive_epsilon_rateH = NumericVector(warmup);
  }
  double sum_epsilon_psi = 0.0;
  NumericVector sum_epsilon_rateP = NumericVector(epsilon_rate_P.size(), 0.0);
  double sum_epsilon_rateH = 0.0;
  
  for (int m = 0; m < M; ++m) {
    
    if (m % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    if (verbose == 1 && m_mod > 0 && m % m_mod == 0) {
      Rcout << "*" ;
    }

    // M-H rate_H
    List out_rate_H = MH_rate_H(data_objects, indolents, prior, 
                                age_at_tau_hp_hats, theta, 
                                epsilon_rate_H, t0);

    theta = out_rate_H["theta"];
    bool accept_rate_H = out_rate_H["accept"];
    double acceptance_prob_rate_H = out_rate_H["probability"];
    
    // M-H rate_P
    List out_rate_P = MH_rate_P(data_objects, indolents, prior, 
                                age_at_tau_hp_hats, theta, 
                                epsilon_rate_P, t0);

    theta = out_rate_P["theta"];
    LogicalVector accept_rate_P = out_rate_P["accept"];
    NumericVector acceptance_prob_rate_P = out_rate_P["probability"];
    
    // Gibbs beta
    theta = gibbs_beta_List(data_objects, prior, age_at_tau_hp_hats, theta, 
                            n_screen_positive_total);

    bool accept_psi;
    double acceptance_prob_psi;
    if (abs(epsilon_psi) > 1e-12) {
      // update (psi, indolent)
      List out_psi_indolent = MH_psi_indolent(data_objects, indolents, prior, 
                                              age_at_tau_hp_hats, theta,
                                              epsilon_psi, t0);

      indolents = out_psi_indolent["indolents"];
      theta = out_psi_indolent["theta"];
      accept_psi = out_psi_indolent["accept"];
      acceptance_prob_psi = out_psi_indolent["probability"];
    } else {
      accept_psi = false;
      acceptance_prob_psi = 0.0;
    }

    // update age_at_tau_hp_hat
    List out_tau  = MH_tau_List(data_objects, indolents, age_at_tau_hp_hats, 
                                theta, t0);

    age_at_tau_hp_hats = out_tau["age_at_tau_hp_hats"];
    List accept_latent = out_tau["accept"];

    // save output
    if (m % thin == 0 && m >= warmup) {
      int start = 0;
      int end;
      IntegerVector ids, itmp;
      LogicalVector ltmp;

      for (int k = 0; k < age_at_tau_hp_hats.size(); ++k) {
        end = start + cases[k] - 1;
        ids = seq(start, end);

        tmp = age_at_tau_hp_hats[k];
        age_at_tau_hp_hat_thin[ids] = tmp;
        
        itmp = indolents[k];
        indolent_thin[ids] = itmp;
        
        ltmp = accept_latent[k];
        accept_latent_thin[ids] = ltmp;
        
        start = end + 1;
      }
      
      ikept += 1;
      
      kRATE_H  [ikept] = theta["rate_H"];
      NumericVector tmp_vec = theta["rate_P"];
      kRATE_P(ikept, _ ) = tmp_vec;
      kPSI     [ikept] = theta["psi"];
      tmp_vec = theta["beta"];
      kBETA(ikept, _ ) = tmp_vec;
      if (save_latent == 1) {
        kage_at_tau_hp_hat  (ikept, _ ) = age_at_tau_hp_hat_thin;
        kINDOLENT(ikept, _ ) = indolent_thin;
        kACCEPT_LATENT(ikept, _ ) = accept_latent_thin;
      } else {
        kage_at_tau_hp_hat  (0, _ ) = age_at_tau_hp_hat_thin;
        kINDOLENT(0, _ ) = indolent_thin;
        kACCEPT_LATENT(0, _ ) = accept_latent_thin;
      }
      kACCEPT_PSI   [ikept] = accept_psi;
      kACCEPT_RATE_H[ikept] = accept_rate_H;
      kACCEPT_RATE_P(ikept, _ ) = accept_rate_P;
    }
    
    if (m < warmup) {
      int shift_to_avg = floor(warmup/2);
      double xi_mp1_wgt = - sqrt(m + 1.0) / gamma / (m + 1.0 + m0);
      double eta = pow(m + 1.0, -kappa);
      
      H_psi += delta - std::min(acceptance_prob_psi, 1.0);
      double xi_mp1 = xi_mp1_wgt * H_psi;
      epsilon_psi = exp((1.0 - eta) * log(epsilon_psi) + eta * xi_mp1);
      adaptive_epsilon_psi[m] = epsilon_psi;
      if (m > shift_to_avg) { sum_epsilon_psi += epsilon_psi; }
      
      for (int i = 0; i < epsilon_rate_P.size(); ++i) {
        H_rateP[i] += delta - std::min(acceptance_prob_rate_P[i], 1.0);
        xi_mp1 = xi_mp1_wgt * H_rateP[i];
        epsilon_rate_P[i] = exp((1.0 - eta) * log(epsilon_rate_P[i]) + eta * xi_mp1);
        if (m > shift_to_avg) { sum_epsilon_rateP[i] += epsilon_rate_P[i]; }
      }
      adaptive_epsilon_rateP(m, _ ) = epsilon_rate_P;
      
      H_rateH += delta - std::min(acceptance_prob_rate_H, 1.0);
      xi_mp1 = xi_mp1_wgt * H_rateH;
      epsilon_rate_H = exp((1.0 - eta) * log(epsilon_rate_H) + eta * xi_mp1);
      adaptive_epsilon_rateH[m] = epsilon_rate_H;
      if (m > shift_to_avg) { sum_epsilon_rateH += epsilon_rate_H; }
      
      if (m == (warmup - 1)) {
        epsilon_psi = sum_epsilon_psi / (m - shift_to_avg);
        epsilon_rate_P = sum_epsilon_rateP / (m - shift_to_avg);
        epsilon_rate_H = sum_epsilon_rateH / (m - shift_to_avg);
      }
    }
    
  } // end runtime
  
  if (verbose == 1 && m_mod > 0) {
    Rcout << std::endl << std::endl;
  }

  // output
  
  List kTHETA = List::create(Named("rate_H") = kRATE_H,
                             Named("rate_P") = kRATE_P,
                             Named("beta") = kBETA,
                             Named("psi") = kPSI);
  List kACCEPT = List::create(Named("rate_H") = kACCEPT_RATE_H,
                              Named("rate_P") = kACCEPT_RATE_P,
                              Named("tau_hp") = kACCEPT_LATENT,
                              Named("psi") = kACCEPT_PSI);
  List epsilon = List::create(Named("rate_H") = epsilon_rate_H,
                              Named("rate_P") = epsilon_rate_P,
                              Named("psi") = epsilon_psi);
  List adaptive_epsilon = List::create(Named("rate_H") = adaptive_epsilon_rateH,
                                       Named("rate_P") = adaptive_epsilon_rateP,
                                       Named("psi") = adaptive_epsilon_psi);
  return List::create(Named("theta") = kTHETA, 
                      Named("tau_hp") = kage_at_tau_hp_hat, 
                      Named("indolent") = kINDOLENT,
                      Named("accept") = kACCEPT,
                      Named("epsilon") = epsilon,
                      Named("adaptive") = adaptive_epsilon,
                      Named("last_theta") = theta); 
}

// [[Rcpp::export]]
NumericVector model_comparison(List data_object, List theta, double t0, int indolent_phase) {

  List prob_tau = compute_prob_tau_obj(data_object, theta, t0);
  NumericVector age_at_tau_hp_hat = rprop_age_at_tau_hp_hat_obj(data_object, prob_tau, theta, t0);
  
  NumericVector prob_indolent(age_at_tau_hp_hat.size());
  IntegerVector indolent(age_at_tau_hp_hat.size());
  if (indolent_phase == 1) {
    prob_indolent = compute_prob_indolent_obj(data_object, theta, age_at_tau_hp_hat);
    indolent = rprop_indolent_obj(data_object, prob_indolent);
  } else {
    prob_indolent.fill(0.0);
    indolent.fill(0);
  }

  NumericVector qlog = dlog_prop_latent_obj(data_object, prob_tau, theta, 
                                            age_at_tau_hp_hat, 
                                            prob_indolent, indolent, t0);

  NumericVector loglik_joint = dlog_likelihood_obj(data_object,
                                                   theta, 
                                                   age_at_tau_hp_hat,
                                                   indolent, t0);
  
  return loglik_joint - qlog;
}
