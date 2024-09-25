#' Bayesian Analysis of Cancer Latency with Auxiliary Variable Augmentation
#'
#' Markov chain Monte Carlo sampler to fit a three-state mixture compartmental
#'  model of cancer natural history to individual-level screening and cancer
#'  diagnosis histories in a Bayesian framework.
#'
#' @details
#'
#' Input \code{theta_0} contains the initial values for all distribution
#'   parameters. The list must include
#'   \itemize{
#'     \item \code{rate_H}: A scalar numeric. The rate for the Weibull
#'     distribution of the healthy compartment.
#'     \item \code{shape_H}: A scalar numeric. The shape parameter for the
#'     Weibull distribution of the healthy compartment.
#'     \item \code{rate_P}: A numeric scalar or named numeric vector. The rate
#'     parameter for each Weibull distribution of the preclinical compartment.
#'     If all participants follow the same Weibull distribution, provide a scalar.
#'     If multiple preclinical Weibull distributions are used, see note below.
#'     \item \code{shape_P}: A scalar numeric. The shape parameter for all
#'     Weibull distributions of the preclinical compartment.
#'     \item \code{beta}: A  scalar numeric or named numeric vector. The assessment
#'     sensitivity. If the sensitivity is the same for all participants, provide
#'     a scalar. If the sensitivity is arm- or screen-type-specific, see note below.
#'     Each element must be in [0, 1].
#'     \item \code{psi}: A scalar numeric. The probability of being indolent.
#'     Must be in [0,1]. If disease is always progressive, this element is
#'     required, but its value must be set to 0.
#'   }
#'
#' Input \code{prior} contains all distribution parameters for the priors.
#'   The list must include
#'   \itemize{
#'     \item \code{rate_P}: A scalar numeric or named vector object. The rate
#'     for the Gamma(shape_P, rate_P) prior on the rate of the Weibull of the
#'     preclinical compartment. If group-specific distributions are used, see
#'     note below.
#'     \item \code{shape_P}: A scalar numeric or named vector object. The shape
#'     for the Gamma(shape_P, rate_P) prior on the rate of the Weibull of the
#'     preclinical compartment. If group-specific distributions are used, see
#'     note below.
#'     \item \code{rate_H}: A scalar numeric. The rate for the
#'     Gamma(shape_H, rate_H) prior on the rate of the Weibull of the healthy
#'     compartment.
#'     \item \code{shape_H}: A scalar numeric. The shape for the
#'     Gamma(shape_H, rate_H) prior on the rate of the Weibull of the healthy
#'     compartment.
#'     \item \code{a_beta}: A positive scalar numeric or named numeric vector.
#'     The first parameter of the Beta(a, b) prior on the assessment sensitivity.
#'     If arm- or screen-type-specific distributions are used, see note below.
#'     If beta is not allowed to change, specify 0.0.
#'     \item \code{b_beta}: A positive scalar numeric or named numeric vector.
#'     The second parameter of the Beta(a, b) prior on the assessment sensitivity.
#'     If arm- or screen-type-specific distributions are used, see note below.
#'     If beta is not allowed to change, specify 0.0.
#'     \item \code{a_psi}: A positive scalar numeric. The first parameter of the
#'     Beta(a, b) prior on the indolence probability. If disease under analysis
#'     does not have an indolent state, this element must be included, but it
#'     will be ignored.
#'     \item \code{b_psi}: A positive scalar numeric. The second parameter of the
#'     Beta(a, b) prior on the indolence probability. If disease under analysis
#'     does not have an indolent state, this element must be included, but it
#'     will be ignored.
#'   }
#'
#' It is possible to assign participants to study arms such that each arm has
#'   its own screening sensitivities and/or rate_P distributions, or to
#'   assign screen-type specific sensitivities.
#'
#'   To designate study arms, each of which will have its own screening
#'   sensitivities:
#'   \itemize{
#'     \item Provide an additional column in \code{data.clinical} named
#'       "arm", which gives the study arm to which each participant is
#'       assigned. For example,
#'       \code{data.clinical$arm = c("Control", "Tx", "Tx", ...)}.
#'     \item Define all beta related prior parameters as named vectors.
#'       For example, \code{prior$a_beta = c("Control" = 1, "Tx" = 38.5)},
#'       and \code{prior$b_beta = c("Control" = 1, "Tx" = 5.8)}
#'     \item Define the initial beta values of theta as a named vector.
#'       For example,  \code{theta_0$beta = c("Control" = 0.75, "Tx" = 0.8)}.
#'   }
#'
#'   Similarly, if using multiple preclinical Weibull distributions
#'     (distributions will have the same shape_P),
#'   \itemize{
#'     \item Provide an additional column in \code{data.clinical} named
#'       "grp.rateP", which assigns each participant to one of the preclinical
#'       Weibull distributions.
#'       For example, \code{data.clinical$grp.rateP = c("rateP1", "rateP2", "rateP2", ... )}.
#'     \item Define the rate_P prior parameter as a named vector.
#'       For example, \code{prior$rate_P <- c("rateP1" = 0.01, "rateP2" = 0.02)}.
#'     \item Define the shape_P prior parameter as a named vector.
#'       For example, \code{prior$shape_P <- c("rateP1" = 1, "rateP2" = 2)}.
#'     \item Define the initial rate_P values of theta as a named vector.
#'       For example,  \code{theta_0$rate_P <- c("rateP1" = 1e-5, "rateP2" = 0.01)}.
#'     \item Define step size of rate_P as a named vector. For example,
#'       \code{epsilon_rate_P <- c("rateP1" = 0.001, "rateP2" = 0.002)}.
#'   }
#'
#'   To assign screen-specific sensitivities,
#'   \itemize{
#'     \item Provide an additional column in \code{data.assess} named
#'       "screen_type", which gives the screening type for each screen. For example,
#'       \code{data.assess$screen_type = c("film", "2D", "2D", ...)}.
#'     \item Define all beta related prior parameters as named vectors.
#'       For example, \code{prior$a_beta = c("film" = 1, "2D" = 38.5)},
#'       and \code{prior$b_beta = c("film" = 1, "2D" = 5.8)}
#'     \item Define the initial beta values of theta as a named vector.
#'       For example,  \code{theta_0$beta = c("film" = 0.75, "2D" = 0.8)}.
#'   }
#'
#'   NOTE: If using integers to indicate group membership, vector names
#'   still must be provided. For example, if group membership is binary 0/1,
#'   vector elements of the prior, initial theta, and step size must be named
#'   as "0" and "1".
#'
#' The adaptive MCMC tuning expression at step m + 1 is defined as
#'   \deqn{\epsilon_{m+1} = (1 - m^{\kappa}) \epsilon_{m} + m^{\kappa}
#'   \xi_{m+1},}{e_{m+1} = (1-m^{kappa}) e_m + m^{kappa} xi_{m+1},}
#'   where
#'   \deqn{\xi_{m+1} = \frac{\sqrt{m}}{\gamma}\frac{1}{m+m_0}
#'   \sum_{i=1}^{m} (\alpha_m - \delta).}{xi_{m+1} = sqrt(m) / gamma *
#'   1/(m+m0) sum(alpha_m - delta).}
#' To initiate the adaptive selection procedure, input \code{adaptive}
#'   must specify the parameters of the above expressions.
#'   Specifically, the provided list must contain elements "delta", the
#'   target acceptance rate; "warmup", the number of iterations to apply step
#'   size correction; and parameters "m0", "kappa", and "gamma".
#'
#' @param data.assess A data.frame. Disease status assessments recorded
#'   during healthy or preclinical compartment, e.g., screenings for disease.
#'   The data must be structured as
#'   \itemize{
#'     \item \code{id}: A character, numeric, or integer object. The unique participant
#'     id to which the record pertains. Multiple records for each id are allowed.
#'     \item \code{age_assess}: A numeric object. The participant's age at time of
#'     assessment.
#'     \item \code{disease_detected}: An integer object. Must be binary 0/1, where
#'     1 indicates that disease was detected at the assessment; 0 otherwise.
#'  }
#'  If the sensitivity parameter (beta) is screen-specific, an additional
#'   column \code{screen_type} is required indicating the type of each
#'   screen.
#' @param data.clinical A data.frame. The clinical data. The data must
#'   be structured as
#'   \itemize{
#'     \item \code{id}: A character, numeric, or integer object. The unique participant
#'       id to which the record pertains. Note these must include those provided in
#'       \code{data.assess}. Must be only 1 record for each participant.
#'     \item \code{age_entry}: A numeric object. The age at time of entry into the study.
#'       Note that this data is used to calculate a normalization; to expedite
#'       numerical integration, it is recommended that the ages be rounded.
#'       Optional input \code{round.age.entry} can be
#'       set to FALSE if this approximation is not desired; however, the
#'       computation time will significantly increase.
#'     \item \code{endpoint_type}: A character object. Must be one of \{"clinical",
#'       "censored", "preclinical"\}. Type "clinical" indicates that disease
#'       was diagnosed in the clinical compartment (i.e., symptomatic). Type
#'       "preclinical" indicates that disease was diagnosed in the preclinical
#'       compartment (i.e., during an assessment). Type
#'       "censored" indicates disease was not diagnosed prior to end of study.
#'     \item \code{age_endpoint}: A numeric object. The participant's age at the
#'       time the endpoint was evaluated.
#'  }
#'  If the sensitivity parameter (beta) is arm-specific, an additional
#'   column \code{arm} is required indicating the study arm to which each
#'   participant is assigned. Similarly, if the preclinical Weibull distribution is
#'   group-specific, an additional column \code{grp.rateP} is required. See Details
#'   for further information.
#' @param baclava.object NULL or a 'baclava' object. To continue a calculation,
#'   provide the object returned by a previous call.
#' @param theta_0 A list object. The initial values for all distribution
#'   parameters. If \code{baclava.object} is a 'baclava' object, this input is ignored.
#'   See Details for further information.
#' @param prior A list object. The prior parameters. If \code{baclava.object} is
#'   a 'baclava' object, this input is ignored. See Details for further
#'   information.
#' @param epsilon_rate_H A small scalar numeric. The Monte Carlo step size
#'   for rate_H (the rate parameter of the Weibull of the healthy compartment).
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param epsilon_rate_P A small scalar numeric or named numeric vector.
#'   The Monte Carlo step size for rate_P (the rate parameter of the Weibull of
#'   the preclinical compartment). If group-specific Weibull distributions
#'   are used, this must be a vector; see Details for further information.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param epsilon_psi A small scalar numeric. The Monte Carlo step size for
#'   parameter psi (the probability of indolence). If disease under
#'   analysis does not have an indolent state, set to 0 and ensure that
#'   the initial value for psi in theta_0 is also 0.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param indolent A logical object. If \code{FALSE}, disease under analysis
#'   does not have an indolent state, i.e., it is always progressive.
#'   This input is provided for convenience; if FALSE, \code{epislon_psi} and
#'   \code{theta_0$psi} will be set to 0.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param t0 A non-negative scalar numeric object. The risk onset age. Must be
#'   less than the earliest assessment age, entry age, and endpoint age.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param M A positive integer object. The number of Monte Carlo samples. This
#'   is the total, i.e., M = adaptive$warmup + n_MCMC.
#' @param thin A positive integer object. Keep each thin-th step of the sampler
#'   after the warmup period, if any, is complete.
#' @param adaptive NULL or named list. If NULL, the step sizes are
#'   not modified in the MCMC. If a list, the parameters for the
#'   adaptive MCMC.
#'   The provided list must contain elements "delta", the target acceptance
#'   rate; "warmup", the number of iterations to apply step size correction;
#'   and parameters "m0", "kappa", and "gamma". See Details for further
#'   information.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param round.age.entry A logical object. If TRUE, the age at time of entry
#'   will be rounded to the nearest integer prior to performing the MCMC.
#'   This data is used to estimate the probability of experiencing clinical
#'   disease prior to entering the study, which is estimated using a
#'   time consuming numerical integration procedure. It is expected that
#'   rounding the ages at time of entry introduces minimal bias. If FALSE,
#'   and ages cannot be grouped, these integrals significantly increase
#'   computation time.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param verbose A logical object. If \code{TRUE}, a progress bar will be shown
#'   during the MCMC.
#' @param save.latent A logical object. If \code{TRUE}, latent variable tau_HP
#'   and indolence will be returned. These can be very large matrices. To
#'   estimate the cohort overdiagnosis probability using \link{cohortODX}(),
#'   this must be set to TRUE.
#'
#' @returns An object of S3 class \code{baclava}, which extends a list object.
#'   \itemize{
#'     \item theta: A list of the posterior distribution parameters at the thinned samples.
#'       \itemize{
#'         \item rate_H: A numeric vector. The rates for the Weibull of the
#'         the healthy compartment.
#'         \item shape_H: A scalar numeric. The input shape_H parameter.
#'         \item rate_P: A numeric matrix. The rates for the Weibull of the
#'          preclinical compartment.
#'         \item shape_P: A scalar numeric. The input shape_P parameter.
#'         \item beta: A numeric matrix. The assessment sensitivities.
#'         \item psi: A numeric vector. The probabilities of indolence. Will be
#'         NA if disease is always progressive.
#'       }
#'     \item tau_hp: If \code{save.latent = TRUE}, a matrix.
#'       The age at time of transition from
#'       healthy to preclinical compartment for each participant at the
#'       thinned samples.
#'     \item indolent: If \code{save.latent = TRUE}, a matrix.
#'       The indolent status for each participant at the
#'       thinned samples. Will be \code{NA} if disease is always progressive.
#'     \item accept: A list of the accept indicator at the thinned samples.
#'       \itemize{
#'         \item rate_H: A numeric vector.
#'         \item rate_P: A numeric matrix.
#'         \item tau_hp: If \code{save.latent = TRUE}, a matrix. Will be NA if
#'           current and new transition ages are Inf.
#'         \item psi: A numeric vector. The probability of indolence. Will be
#'         \code{NA} if disease is always progressive.
#'       }
#'     \item epsilon: A list. The step sizes for each parameter.
#'     \item adaptive: A list. Settings for the adaptive procedure. Will be NA
#'       if adaptive procedure not requested.
#'     \item last_theta: A list. The theta parameters of the last MCMC iteration.
#'     \item prior: A list. The provided parameters of the prior distributions.
#'     \item setup: A list of inputs provided to the call.
#'       \itemize{
#'         \item t0: The input age of risk onset.
#'         \item indolent: TRUE if disease is not progressive.
#'         \item round.age.entry: TRUE if age at entry was rounded to the nearest
#'           whole number.
#'         \item groups.beta: A vector of the beta grouping values.
#'         \item groups.rateP: A vector of the rate_P grouping values.
#'         \item thin: The number of samples dropped between kept MCMC iterations.
#'         \item initial.theta: theta_0 as provided by user.
#'         \item initial.prior: prior as provided by user.
#'      }
#'     \item clinical.groupings: A data.frame of the original data's arm/rateP grouping.
#'     \item screen_types: A data.frame of the original data's screen type grouping.
#'     \item call: The matched call.
#'   }
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
#' # This is for illustration only -- the number of Gibbs samples should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        theta_0 = theta_0,
#'                        prior = prior)
#'
#' summary(example)
#' print(example)
#'
#' # To continue this calculation
#' example_continued <- fit_baclava(data.assess = data.screen,
#'                                  data.clinical = data.clinical,
#'                                  baclava.object = example)
#' @importFrom dplyr count
#' @importFrom tidyr nest
#' @import Rcpp
#' @import RcppNumerical
#' @include computeEndpoints.R RcppExports.R utilities.R
#' @useDynLib baclava
#' @export
fit_baclava <- function(data.assess,
                        data.clinical,
                        baclava.object = NULL,
                        M = 100L,
                        thin = 1L,
                        t0 = 0,
                        theta_0 = list(),
                        prior = list(),
                        epsilon_rate_H = 1e-3,
                        epsilon_rate_P = 1e-3,
                        epsilon_psi = 1e-3,
                        indolent = TRUE,
                        adaptive = NULL,
                        round.age.entry = TRUE,
                        verbose = TRUE,
                        save.latent = FALSE) {

  if (!is.null(baclava.object)) {
    # this is a continuation of a previous call
    # extract all model specification information from previous object
    # list returned sets t0, indolent, theta_0, prior, epsilon_rate_H,
    #   epsilon_rate_P, and epsilon_psi
    tmp_list <- .continuation(baclava.object, verbose = TRUE)
    for (i in seq_along(tmp_list)) assign(names(tmp_list)[i], tmp_list[[i]])
    # adaptive procedure is not allowed in a continuation
    adaptive <- NULL
  }

  stopifnot(
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`theta_0` must be a list containing rate_H, shape_H, rate_P, shape_P, beta, and psi" =
      is.vector(theta_0, mode = "list") &&
      all(c("rate_H", "shape_H", "rate_P", "shape_P", "beta", "psi") %in% names(theta_0)),
    "`prior` must be a list containing rate_H, shape_H, rate_P, shape_P, a_beta, b_beta, a_psi, and b_psi" =
      is.list(prior) &&
      all(c("rate_H", "shape_H", "rate_P", "shape_P", "a_beta", "b_beta", "a_psi", "b_psi") %in% names(prior)),
    "`epsilon_rate_H` must be a positive scalar" =
      is.vector(epsilon_rate_H, mode = "numeric") &&
      length(epsilon_rate_H) == 1L && epsilon_rate_H >= 0.0,
    "`epsilon_rate_P` must be a positive scalar or vector of positive values" =
      is.vector(epsilon_rate_P, mode = "numeric") && all(epsilon_rate_P >= 0.0),
    "`epsilon_psi` must be a scalar numeric" = is.vector(epsilon_psi, mode = "numeric") &&
      length(epsilon_psi) == 1L && epsilon_psi >= 0.0,
    "`indolent` must be logical" = is.vector(indolent, mode = "logical") &&
      length(indolent) == 1L,
    "`t0` must be a non-negative scalar numeric" =
      is.vector(t0, mode = "numeric") && length(t0) == 1L && t0 >= 0.0,
    "`M` must be a positive integer" = is.vector(M, mode = "numeric") &&
      length(M) == 1L && isTRUE(all.equal(M, as.integer(M))) && M > 0.0,
    "`thin` must be a positive integer" = is.vector(thin, mode = "numeric") &&
      length(thin) == 1L && isTRUE(all.equal(thin, as.integer(thin))) && thin > 0.0,
    "`adaptive` must be NULL or a list" = is.null(adaptive) || is.list(adaptive),
    "`round.age.entry` must be logical" = is.logical(round.age.entry),
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") &&
      length(verbose) == 1L,
    "`baclava.object` must be NULL or of class 'baclava'" = is.null(baclava.object) ||
      "baclava" %in% class(baclava.object)
  )

  provided_theta <- theta_0
  provided_prior <- prior

  ### Verify data.frame inputs

  clinical.ids <- as.character(data.clinical$id)
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
                                   prior = prior,
                                   epsilon = epsilon_rate_P)

  data.clinical$grp.rateP <- tmp_list$grp.rateP
  prior <- tmp_list$prior
  theta_0 <- tmp_list$theta
  epsilon_rate_P <- tmp_list$epsilon
  groups_rateP <- levels(data.clinical$grp.rateP)
  data.clinical$grp.rateP <- unclass(data.clinical$grp.rateP)

  ### Verify adaptive input

  if (is.null(adaptive)) {
    adaptive = list("warmup" = 0L,
                    "delta" = 0.0,
                    "gamma" = 0.0,
                    "kappa" = 0.0,
                    "m0" = 0L)
  } else {
    .testAdaptiveInput(adaptive)
    if (adaptive$warmup >= M) {
      stop("M cannot be <= adaptive$warmup", call. = FALSE)
    }
  }

  adaptive$warmup <- as.integer(adaptive$warmup)
  adaptive$m0 <- as.integer(adaptive$m0)


  ### Verify theta_0
  theta_0 <- .testThetaInput(theta_0)

  ### Verify prior
  prior <- .testPriorInput(prior)

  ### misc data prepping

  n_participants <- nrow(data.clinical)

  # for each screen-type, number of screens that detected disease
  n_screen_positive_total <-
    table(data.assess$disease_detected,
          data.assess$screen_type)["1", ] |> unname()

  if (round.age.entry) {
    message("NOTE: age at entry rounded to nearest whole number")
    data.clinical$age_entry <- as.integer(round(data.clinical$age_entry, 0L))
  }

  data.objs <- .makeDataObjectsFull(data.clinical = data.clinical,
                                    data.assess = data.assess,
                                    groups.rateP = groups_rateP,
                                    t0 = t0)

  # it is convenient to have the rateP vector stored outside of the individual
  # data objects
  theta_0$irateP <- lapply(data.objs, "[[", "irateP") |> unlist()

  # returned data will be reordered to match the order of data as provided
  orig_order <- match(
    lapply(data.objs, "[[", "ids") |> unlist(),
    clinical.ids
  )
  if (any(is.na(orig_order))) {
    warning("problem matching ids for reordering; order of results may not match order of provided data.clinical", call. = FALSE)
    orig_order <- seq_len(nrow(data.clinical))
  }

  # initialization of latent variables -- these functions are defined in C++
  theta_0$scale_H <- theta_0$rate_H^{-1.0 / theta_0$shape_H}
  theta_0$scale_P <- theta_0$rate_P^{-1.0 / theta_0$shape_P}

  prob_tau_0 <- compute_prob_tau_List(data.objs, theta_0, t0)

  if (!is.null(baclava.object)) {
    # extract last age_at_tau_hp and indolents from prior object
    tmp_list <- .getPreviousResults(data.objects = data.objs,
                                    baclava.object = baclava.object,
                                    indolent = indolent)
    age_at_tau_hp_hats = tmp_list$age.at.tau.hp.hats
    indolents = tmp_list$indolents
  } else {
    # calculate initial values for age_at_tau_hp and indolents
    age_at_tau_hp_hats <- rprop_age_at_tau_hp_hat_List(data.objs, prob_tau_0, theta_0, t0)
    prob_indolent_0 <- compute_prob_indolent_List(data.objs, age_at_tau_hp_hats, theta_0)
    indolents <- rprop_indolent_List(data.objs, prob_indolent_0)
  }

  if (!indolent || {epsilon_psi < 1e-12 && theta_0$psi < 1e-12}) {
    # if disease is always progressive, adjust parameters and set indolents to
    # 0 for all cases
    indolent <- FALSE
    epsilon_psi <- 0.0
    theta_0$psi <- 0.0

    indolents <- list()
    for (i in seq_along(data.objs)) {
      indolents[[i]] <- rep(0L, data.objs[[i]]$n)
    }
  }

  out_object <- MCMC_cpp_internal(data_objects = data.objs,
                                  indolents = indolents,
                                  prior = prior,
                                  age_at_tau_hp_hats = age_at_tau_hp_hats,
                                  theta = theta_0,
                                  epsilon_rate_H = epsilon_rate_H,
                                  epsilon_rate_P = epsilon_rate_P,
                                  epsilon_psi = epsilon_psi,
                                  t0 = t0,
                                  M = M,
                                  thin = thin,
                                  M_thin = (M - max(adaptive$warmup, 0)) / thin,
                                  n_obs = n_participants,
                                  n_screen_positive_total = n_screen_positive_total,
                                  adaptive = adaptive,
                                  verbose = as.integer(verbose),
                                  save_latent = as.integer(save.latent))

  # re-order the columns of output, so they match the order of the input
  out_object$tau_hp[, orig_order] = out_object$tau_hp
  out_object$indolent[, orig_order] = out_object$indolent
  out_object$accept$tau_hp[, orig_order] = out_object$accept$tau_hp

  # These need to be named to ensure continuation is in correct order
  colnames(out_object$tau_hp) <- clinical.ids
  colnames(out_object$indolent) <- clinical.ids
  colnames(out_object$accept$tau_hp) <- clinical.ids

  out_object$theta$shape_H <- theta_0$shape_H
  out_object$theta$shape_P <- theta_0$shape_P

  colnames(out_object$theta$rate_P) <- groups_rateP
  colnames(out_object$theta$beta) <- screen_types
  names(out_object$last_theta$rate_P) <- groups_rateP
  names(out_object$last_theta$beta) <- screen_types

  out_object$prior <- prior
  out_object$setup <- list("t0" = t0,
                           "indolent" = indolent,
                           "round.age.entry" = round.age.entry,
                           "groups.beta" = screen_types,
                           "groups.rateP" = groups_rateP,
                           "thin" = thin,
                           "initial.theta" = provided_theta,
                           "initial.prior" = provided_prior)

  out_object$groupings <- data.clinical[, c("id", "arm", "grp.rateP")]
  data.assess$screen.id <- sapply(table(data.assess$id), seq, from = 1L, by = 1L) |> unlist()
  out_object$screen_types <- data.assess[, c("id", "screen.id", "screen_type")]

  if (!indolent) {
    out_object$indolent <- NA_real_
    out_object$accept$psi <- NA_real_
    out_object$theta_0$psi <- NA_real_
  }

  out_object$call <- match.call()

  class(out_object) <- "baclava"

  out_object
}
