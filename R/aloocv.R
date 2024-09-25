#' Approximate Leave-One-Out Cross-Validation
#'
#' Approximate leave-one-out cross-validation computed from the posterior
#'   draws of the Markov chain Monte Carlo sampler as implemented in
#'   \link{fit_baclava}().
#'
#' Computes the predictive fit of a model.
#'   For each individual and each MCMC draw, the function approximates the
#'   marginal likelihood via importance sampling. It samples J.increment values
#'   of the individual's latent variables using the Metropolis-Hastings proposal
#'   distributions and computes the effective sample size (ESS) of the
#'   importance sampling procedure. If the target ESS is not met, J.increment
#'   additional samples are taken, and the ESS is re-evaluated.
#'   This is repeated until either the ESS is
#'   satisfied or J.max samples have been drawn.
#'
#' @param object The value object returned by \link{fit_baclava}().
#' @param data.clinical A data.frame object. The clinical data on which the
#'   model is assessed. The data must be structured as for \code{fit_baclava()};
#'   specifically, it must contain
#'   \itemize{
#'     \item{id:} A character, numeric, or integer object. The unique participant
#'       id to which the record pertains. Note these must include those provided in
#'       \code{data.assess}. Must be only 1 record for each participant.
#'     \item{age_entry:} A numeric object. The age at time of entry into the study.
#'       Note that this data is used to calculate the normalization; to expedite
#'       numerical integration, it is recommended that the ages be rounded to
#'       minimize repeated calculations. Optional input `round.age.entry` can be
#'       set to FALSE if this approximation is not desired; however, the
#'       computation time will significantly increase.
#'     \item{endpoint_type:} A character object. Must be one of \{"clinical",
#'       "censored", "preclinical"\}. Type "clinical" indicates that disease
#'       was diagnosed in the clinical compartment (i.e., symptomatic). Type
#'       "preclinical" indicates that disease was diagnosed in the pre-clinical
#'       compartment (i.e., during an assessment). Type
#'       "censored" indicates participant was censored.
#'     \item{age_endpoint:} A numeric object. The participant's age at the
#'       time the endpoint was evaluated.
#'  }
#'  If the sensitivity parameter (beta) is arm-specific, an additional
#'   column \code{arm} is required indicating the study arm to which each
#'   participant is assigned. Similarly, if the preclinical Weibull distribution is
#'   group-specific, an additional column \code{grp.rateP} is required. See Details
#'   for further information.
#'
#' @param data.assess A data.frame object. The disease status assessment data
#'   on which the model is assessed. The data must be structured as for
#'   \code{fit_baclava()}; specifically, the data must contain
#'   \itemize{
#'     \item{id:} A character, numeric, or integer object. The unique participant
#'     id to which the record pertains.
#'     \item{age_assess:} A numeric object. The participant's age at time of
#'     assessment.
#'     \item{disease_detected:} An integer object. Must be binary 0/1, where
#'     1 indicates that disease was detected at the assessment; 0 otherwise.
#'  }
#'  If the sensitivity parameter (beta) is screen-type specific, an additional
#'   column \code{screen_type} is required indicating the type of each
#'   screen.
#'
#' @param J.increment An integer object. The number of replicates of each
#'   participant to generate in each iteration of the importance sampling procedure
#'   to attain desired effective sample size.
#' @param J.max An integer object. The maximum number of samples to be
#'   drawn.
#' @param ess.target An integer object. The target effective sample size in the
#'   importance sampling procedure.
#' @param n.core An integer object. The function allows for the outer loop
#'   across participants to be run in parallel using foreach().
#' @param verbose A logical object. If TRUE, progress information will be
#'   printed. This input will be ignored if n.core > 1.
#' @param lib An optional character vector allowing for library path
#'   to be provided to cluster.
#'
#' @returns A list object. Element \code{summary} contains the min, mean, and
#'   the 1% quantile of the ESS; the max, mean, and 99% quantile of J; the
#'   likelihood; and the individual-level and estimated predictive fit.
#'   Element \code{result} contains the likelihood, ESS, and J for each
#'   MCMC sample for each participant.
#'
#' @examples
#'
#' data(screen_data)
#'
#' theta_0 <- list("rate_H" = 7e-4, "shape_H" = 2.0,
#'                 "rate_P" = 0.5  , "shape_P" = 1.0,
#'                 "beta" = 0.9, psi = 0.4)
#' prior <- list("rate_H" = 0.01, "shape_H" = 1,
#'               "rate_P" = 0.01, "shape_P" = 1,
#'               "a_psi" = 1/2 , "b_psi" = 1/2,
#'               "a_beta" = 38.5, "b_beta" = 5.8)
#'
#' # This is for illustration only -- the number of MCMC samples should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        theta_0 = theta_0,
#'                        prior = prior,
#'                        thin = 10L)
#'
#' res <- aloocv(example, data.clinical, data.screen)
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
aloocv <- function(object,
                   data.clinical, data.assess,
                   J.increment = 75L,
                   J.max = 225L,
                   ess.target = 50L,
                   n.core = 1L,
                   verbose = TRUE,
                   lib = NULL) {

  stopifnot(
    "`object` must be of class 'baclava'" = !missing(object) &&
      "baclava" %in% class(object),
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`J.increment` must be a positive integer" = is.vector(J.increment, mode = "numeric") &&
      length(J.increment) == 1L && isTRUE(all.equal(J.increment, round(J.increment))) &&
      J.increment > 0,
    "`J.max` must be a positive integer >= J.increment" = is.vector(J.max, mode = "numeric") &&
      length(J.max) == 1L && isTRUE(all.equal(J.max, round(J.max))) &&
      J.max >= J.increment,
    "`ess.target` must be a positive integer" = is.vector(ess.target, mode = "numeric") &&
      length(ess.target) == 1L && isTRUE(all.equal(ess.target, round(ess.target))) &&
      ess.target > 0,
    "`n.core` must be a positive integer" = is.vector(n.core, mode = "numeric") &&
      length(n.core) == 1L && isTRUE(all.equal(n.core, round(n.core))) &&
      n.core > 0,
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") &&
      length(verbose) == 1L,
    "`lib` must be a character vector" = is.null(lib) || is.character(lib)
  )

  t0 <- object$setup$t0
  shape_H <- object$theta$shape_H
  shape_P <- object$theta$shape_P

  theta_0 <- object$setup$initial.theta
  prior <- object$setup$initial.prior

  n_mcmc_samples <- length(object$theta$rate_H)
  n_participants <- nrow(data.clinical)

  data.assess <- .testDataAssess(data.assess)
  data.clinical <- .testDataClinical(data.clinical)

  # all ids in data.assess must be in data.clinical
  # note this is not true in the other direction
  if (!all(data.assess$id %in% data.clinical$id)) {
    stop("participant ids in `data.assess` not found in `data.clinical`",
         call. = FALSE)
  }

  ### Arm- or Screen-specific beta

  tmp_list <- .processSensitivity(data.assess, data.clinical, theta_0, prior)
  data.assess$screen_type <- tmp_list$screen.type
  prior <- tmp_list$prior
  theta_0 <- tmp_list$theta
  data.clinical$arm <- tmp_list$arm

  ### Verify Screen-Type inputs

  # this procedure converts screen_type to a factor and ensures that
  # theta and prior have appropriately named elements
  # theta and prior are returned in the ordering of the factor levels,
  # thus unclassing the factor variable gives the index of prior and
  # theta pertaining to that specific screen type
  tmp_list <- .testGroupBetaInput(screen.type = data.assess$screen_type,
                                  theta = theta_0,
                                  prior = prior)

  data.assess$screen_type <- tmp_list$screen.type
  prior <- tmp_list$prior
  theta_0 <- tmp_list$theta
  screen_types <- levels(data.assess$screen_type)
  data.assess$screen_type <- unclass(data.assess$screen_type)

  ### Verify rateP inputs

  # this procedure converts rate_p to a factor and ensures that
  # theta and prior have appropriately named elements
  # theta and prior are returned in the ordering of the factor levels,
  # thus unclassing the factor variable gives the index of prior and
  # theta pertaining to that specific rateP group
  if ("grp.rateP" %in% colnames(data.clinical)) {
    tmp <- data.clinical$grp.rateP
  } else {
    tmp <- NULL
  }
  tmp_list <- .testGroupRatePInput(grp.rateP = tmp,
                                   n = nrow(data.clinical),
                                   theta = theta_0,
                                   prior = prior)

  data.clinical$grp.rateP <- tmp_list$grp.rateP
  prior <- tmp_list$prior
  theta_0 <- tmp_list$theta
  groups_rateP <- levels(data.clinical$grp.rateP)
  data.clinical$grp.rateP <- unclass(data.clinical$grp.rateP)

  if (object$setup$round.age.entry) {
    data.clinical$age_entry <- as.integer(round(data.clinical$age_entry, 0L))
  }

  data.objs <- .makeDataObjectsFull(data.clinical = data.clinical,
                                    data.assess = data.assess,
                                    groups.rateP = groups_rateP,
                                    t0 = t0)

  theta_0$irateP <- lapply(data.objs, "[[", "irateP") |> unlist()

  if (n.core > 1L) {
    verbose <- FALSE
    cl <- makeCluster(n.core)
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    if (!is.null(lib)) clusterCall(cl, function() { .libPaths(lib) })
    clusterEvalQ(cl, library(doParallel))
  }

  i <- NULL # to circumvent undefined global variable warnings

  foreach_result <- foreach(i = seq_along(data.clinical$id)) %dopar% {

    if (verbose && {{i %% floor(nrow(data.clinical)*0.05)} < 0.5}) {
      print(c(i, nrow(data.clinical)))
    }

    J <- integer(n_mcmc_samples)
    ess <- numeric(n_mcmc_samples)
    lik <- numeric(n_mcmc_samples)
    convergence_flag <- logical(n_mcmc_samples)

    # data of i-th subject
    data.clinical_i = data.clinical[i, ]
    data.assess_i   = data.assess[data.assess$id %in% data.clinical_i$id, ]

    # replicate the data J_increment times
    data.clinical_i_rep = data.frame("id" = 1:J.increment,
                                     "age_entry" = rep(data.clinical_i$age_entry, J.increment),
                                     "endpoint_type" = rep(data.clinical_i$endpoint_type, J.increment),
                                     "age_endpoint" = rep(data.clinical_i$age_endpoint, J.increment),
                                     "arm" = rep(data.clinical_i$arm, J.increment),
                                     "grp.rateP" = rep(data.clinical_i$grp.rateP, J.increment))

    data.assess_i_rep = data.frame("id" = rep(1:J.increment, each=nrow(data.assess_i)),
                                   "age_assess" = rep(data.assess_i$age_assess, J.increment),
                                   "disease_detected" = rep(data.assess_i$disease_detected, J.increment),
                                   "screen_type" = rep(data.assess_i$screen_type, J.increment))

    # extract ages for each participant as a list
    data.objs <- suppressMessages(.makeDataObjectsFull(data.clinical_i_rep,
                                                       data.assess_i_rep,
                                                       groups.rateP = groups_rateP,
                                                       t0, all.types = FALSE)[[1L]])

    for (s in 1L:n_mcmc_samples) {

      theta_sample = list("rate_H" = object$theta$rate_H[s],
                          "shape_H" = shape_H,
                          "rate_P" = object$theta$rate_P[s, ],
                          "shape_P" = shape_P,
                          "beta" = object$theta$beta[s, ],
                          "psi" = object$theta$psi[s])
      theta_sample$irateP <- theta_0$irateP
      theta_sample$scale_H <- theta_sample$rate_H^{-1.0 / theta_sample$shape_H}
      theta_sample$scale_P <- theta_sample$rate_P^{-1.0 / theta_sample$shape_P}
      if (!object$setup$indolent) theta_sample$psi <- 0.0

      integrand <- NULL
      n_increment <- 0L

      while(any(ess[s] < ess.target) && n_increment * J.increment < J.max) {

        integrand_cpp <- model_comparison(data.objs, theta_sample, t0,
                                          as.integer(object$setup$indolent))

        # IS estimate
        integrand <- c(integrand, exp(unlist(integrand_cpp)))
        lik[s] <- mean(integrand)

        # ESS
        integrand_norm <- integrand / sum(integrand)
        ess[s] <- 1.0 / sum(integrand_norm^2)

        n_increment <- n_increment + 1L
      }
      J[s] = n_increment * J.increment
    }

    list("i" = i, "lik" = lik, "ess" = ess, "J" = J)

  }

  i_ordering <- lapply(foreach_result, "[[", "i") |> unlist()
  lik_all <- lapply(foreach_result, "[[", "lik") |> do.call(what = rbind)
  ess_all <- lapply(foreach_result, "[[", "ess") |> do.call(what = rbind)
  J_all <- lapply(foreach_result, "[[", "J") |> do.call(what = rbind)

  if (any(J_all == J.max & ess_all < ess.target)) {
    warning("Target ESS not met for ",
            sum(J_all == J.max & ess_all < ess.target),
            " MCMC sample(s)/participants",
            call. = FALSE)
  }
  results <- list(shape_H = shape_H,
                  shape_P = shape_P,
                  lik = lik_all[i_ordering, ],
                  ess = ess_all[i_ordering, ],
                  J = J_all[i_ordering, ])

  results_summary <- list("ess_min" = min(ess_all),
                          "ess_mean" = mean(ess_all),
                          "ess_q001" = quantile(ess_all, probs = 0.01),
                          "J_max" = max(J_all),
                          "J_mean" = mean(J_all),
                          "J_q099" = quantile(J_all, probs = 0.99),
                          "lpd" = rowMeans(results$lik) |> log() |> sum(),
                          "elpd_loo_vec" = log(1.0 / rowMeans(1.0 / results$lik)))
  results_summary$elpd_loo <- mean(results_summary$elpd_loo_vec)

  list("result" = results, "summary" = results_summary)
}