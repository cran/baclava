#' Probability of death in the interval [t0, t1]
#'
#' @noRd
#' @param t0 A numeric scalar. The lower boundary of the time interval in units
#'   equivalent to those used to define \code{prob}.
#' @param t1 A numeric scalar. The upper boundary of the time interval in units
#'   equivalent to those used to define \code{prob}.
#' @param otherCauses A data.frame object. Must contain 2 columns with headers
#'  Age, the age at which the corresponding probability pertains and
#'  Rate, the probability of death from other causes.
#'
#' @returns A scalar numeric. The probability of death from other causes in the
#'   specified interval.
#'   
#' @keywords internal
computeProbabilityDeath <- function(t0, t1, otherCauses) {

  interval_containing_min <- max(findInterval(t0, otherCauses$Age), 1L)
  interval_containing_max <- findInterval(t1, otherCauses$Age, left.open = TRUE)

  if (interval_containing_min != interval_containing_max) {
    intervals_spanned <- interval_containing_min:interval_containing_max
    ages <- c(otherCauses$Age[intervals_spanned], t1) |> unique()
    ages[1L] <- max(t0, ages[1L])
  } else {
    ages <- c(t0, t1)
    intervals_spanned <- interval_containing_min
  }

  # vector of probabilities of dying in each interval
  prob_death_interval <- diff(ages) * otherCauses$Rate[intervals_spanned]

  # vector of probabilities of surviving in each interval
  prob_survival_interval <- 1.0 - prob_death_interval

  # probability of surviving all the intervals
  prob_survival <- prod(prob_survival_interval)

  # probability of dying in any interval
  1.0 - prob_survival
}

#' Cohort-specific Overdiagnosis
#'
#' Estimates the overall and screening specific overdiagnosis probability
#'   for the cohort of the original analysis.
#'
#' @param object A 'baclava' object. The value object returned by \link{fit_baclava}().
#' @param data.assess A data.frame object. Disease status assessments recorded
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
#'  This input should be identical to that provided to obtain \code{object}.
#' @param data.clinical A data.frame object. The clinical data. The data must
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
#'  This input should be identical to that provided to obtain \code{object}.
#' @param other.cause.rates A data.frame object. Age specific incidence rates
#'   that do not include the disease of interest. Must contain columns "Rate"
#'   and "Age".
#' @param plot A logical object. If TRUE, generates a boxplot of the overdiagnosis
#'   probability for each individual as a function of the screen at which
#'   disease was detected. Includes only the consecutive screens for which more
#'   than 1\% of the screen detected cases were detected.
#'
#' @returns A list object. 
#' \itemize{
#'   \item \code{all} An n x S matrix containing the estimated overdiagnosis
#'     probability for each individual (n) and each posterior parameter set (S).
#'   \item \code{mean.individual} A vector containing the mean across S of
#'     the estimated overdiagnosis for each individual, i.e., \code{rowMeans(all)}.
#'   \item \code{mean.overall} A numeric, the mean overdiagnosis
#'     probability across all posterior parameter sets and screen-detected cases, 
#'     i.e., \code{mean(all)}. 
#'   \item \code{summary.by.screen} A matrix containing the summary statistics of 
#'     \code{mean.individual} for the individuals detected positive at
#'     each screen, i.e., \code{summary(mean.individual[diagnosis_screen_id == i])}. 
#' }
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
#'                        prior = prior,
#'                        save.latent = TRUE)
#'
#' # if rates are not available, an all cause dataset is provided in the package
#' # NOTE: these predictions will be over-estimated
#'
#' data(all_cause_rates)
#' all_cause_rates <- all_cause_rates[, c("Age", "both")]
#' colnames(all_cause_rates) <- c("Age", "Rate")
#'
#' cohort_odx <- cohortODX(object = example,
#'                         data.clinical = data.clinical,
#'                         data.assess = data.screen,
#'                         other.cause.rates = all_cause_rates,
#'                         plot = FALSE)
#'
#' @importFrom stats rweibull
#' @include utilities.R
#' @export
cohortODX <- function(object, data.clinical, data.assess,
                      other.cause.rates = NULL,
                      plot = TRUE) {

  stopifnot(
    "`object` must be an object of class 'baclava'" = !missing(object) &&
      inherits(object, "baclava"),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`other.cause.rates` must be NULL or a data.frame" =
      {is.null(other.cause.rates) || {is.data.frame(other.cause.rates) &&
          all(c("Rate", "Age") %in% colnames(other.cause.rates))}},
    "`plot` must be logical" = is.vector(plot, mode = "logical") && length(plot) == 1L
  )
  
  if (nrow(object$tau_hp) == 1L) {
    stop("cannot calculate cohortODX -- latent data not saved at each iteration\n\t",
         "fit_baclava analysis must be performed with save.latent = TRUE",
         call. = FALSE)
  }
    
  # if the same data as provided in original analysis, ordering will be the same
  # after this
  data.assess <- .testDataAssess(data.assess)
  data.clinical <- .testDataClinical(data.clinical)

  if (nrow(data.clinical) != nrow(object$groupings)) {
    stop('data.clinical does not have the same dimensions as the analysis cohort',
         call. = FALSE)
  }

  if (nrow(data.assess) != nrow(object$screen_types)) {
    stop('data.assess does not have the same dimensions as the analysis cohort',
         call. = FALSE)
  }

  # because of testing call above -- these should be in the correct order, but
  # just in case!
  idx <- match(object$groupings$id, data.clinical$id)
  if (any(is.na(idx))) {
    stop("unable to align analysis object with provided data; verify ids")
  }
  data.clinical$grp.rateP <- object$groupings$grp.rateP[idx]

  # Take the subset of screen-detected individuals
  screen_detected <- data.clinical$endpoint_type == "preclinical"
  n_screen_detected <- sum(screen_detected)

  data.clinical <- data.clinical[screen_detected, ]
  data.assess <- data.assess[data.assess$id %in% data.clinical$id, ]

  n_screens <- table(data.assess$id)

  scale_p_matrix <- object$theta$rate_P^(-1.0 / object$theta$shape_P)
  
  idx <- match(data.clinical$id, colnames(object$tau_hp))
  if (any(is.na(idx))) {
    stop("unable to align column names of baclava results with data ids", call. = FALSE)
  }
  
  Z_tau_screen <- object$tau_hp[, idx, drop = FALSE]
  Z_I_screen <- object$indolent[, idx, drop = FALSE]
  
  # containers
  S <- nrow(scale_p_matrix)
  prob_odx <- matrix(0.0, nrow = n_screen_detected, ncol = S)

  for (i in seq_len(n_screen_detected)) {
    irateP <- data.clinical$grp.rateP[i]
    for (s in seq_len(S)) {
      if (Z_I_screen[s, i] == 1L) {
        prob_odx[i, s] <- 1.0
        next
      }

      tau_PC <- 0.0
      while (tau_PC < data.clinical$age_endpoint[i]) {
        tau_PC <- Z_tau_screen[s, i] +
          stats::rweibull(1, object$theta$shape_P, scale_p_matrix[s, irateP])
      }

      prob_odx[i, s] <- computeProbabilityDeath(t0 = data.clinical$age_endpoint[i],
                                                t1 = tau_PC,
                                                other.cause.rates)
    }
  }
  
  rownames(prob_odx) <- data.clinical$id

  # overall ODX rate
  prob_odx <- prob_odx * 100.0
  overall <- sum(prob_odx) / {S * n_screen_detected}

  prob_by_individual <- rowMeans(prob_odx)

  by_screen <- sapply(1L:max(n_screens),
                      function(i) summary(prob_by_individual[n_screens == i]))
  by_screen <- data.frame(t(by_screen))
  by_screen$n <- tabulate(n_screens)

  if (plot) {
    five_percent <- ceiling(n_screen_detected * 0.01)
    large_enough <- which(tabulate(n_screens) >= five_percent)
    large_enough <- large_enough[c(0L, diff(large_enough)) <= 1L]

    df <- data.frame(prob_odx = prob_by_individual,
                     screen = factor(n_screens))

    gg <- ggplot(subset(df, df$screen %in% large_enough),
                 aes(.data$screen, .data$prob_odx)) +
      geom_boxplot() +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "errorbar", color = "red", width = 0.75,
                   show.legend = TRUE) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)) +
      ylab("Overdiagnosis Probability") +
      xlab("Screen at which Disease Detected") +
      ggtitle("Overdiagnosis Probability vs Screen at which Disease Detected")
    
    print(gg)
  }

  list("all" = prob_odx, 
       "mean.individual" = prob_by_individual,
       "mean.overall" = overall, 
       "summary.by.screen" = by_screen)
  
}