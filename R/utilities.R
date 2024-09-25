.continuation <- function(baclava.object, verbose) {
  
  if (verbose) message("\n", rep("* ", 10), "\n")
  if (verbose) message("Continuation of previous run")
  
  lst <- list()
  
  lst$t0 <- baclava.object$setup$t0
  lst$theta_0 <- baclava.object$last_theta
  lst$prior <- baclava.object$prior
  lst$epsilon_rate_H <- baclava.object$epsilon$rate_H
  lst$epsilon_rate_P <- baclava.object$epsilon$rate_P
  lst$round.age.entry <- baclava.object$setup$round.age.entry
  lst$thin <- baclava.object$setup$thin
  lst$indolent <- baclava.object$setup$indolent
  
  if (lst$indolent) {
    lst$epsilon_psi <- baclava.object$epsilon$psi
  } else {
    lst$epsilon_psi <- 0.0
    lst$theta_0$psi <- 0.0
  }
  
  if (verbose) {
    message("\t t0 = ", lst$t0, "\t indolent = ", lst$indolent)
    message("\t theta_0 = list(rate_H = ", lst$theta_0$rate_H, ", shape_H = ", lst$theta_0$shape_H, ", \n",
            "\t                rate_P = ", lst$theta_0$rate_P, ", shape_P = ", lst$theta_0$shape_P, ", \n",
            "\t                beta = ", lst$theta_0$beta, ", \n",
            "\t                psi = ", lst$theta_0$psi, ")")
    message("\t prior = list(rate_H = ", lst$prior$rate_H, ", shape_H = ", lst$prior$shape_H, ", \n",
            "\t              rate_P = ", lst$prior$rate_P, ", shape_P = ", lst$prior$shape_P, ", \n",
            "\t              a_beta = ", lst$prior$a_beta, ", b_beta = ", lst$prior$b_beta, ", \n",
            "\t              a_psi = ", lst$prior$a_psi, ", b_psi = ", lst$prior$b_psi, ")")
    message("\t epsilon_rate_H = ", lst$epsilon_rate_H)
    message("\t epsilon_rate_P = ", lst$epsilon_rate_P)
    message("\t epsilon_psi = ", lst$epsilon_psi)
    message("\n", rep("* ", 10), "\n")
  }
  lst
}

.processSensitivity <- function(data.assess, data.clinical, theta, prior) {
  
  if (is.vector(theta$beta, "numeric")) {
    
    if (length(theta$beta) == 1L) {
      # if beta is a scalar numeric, all screens have the same sensitivity
      
      if (!{"screen_type" %in% colnames(data.assess)}) {
        # if no screen-type designation provided, set to a default value
        screen_type <- rep("screen", nrow(data.assess))
        names(theta$beta) <- "screen"
        names(prior$a_beta) <- "screen"
        names(prior$b_beta) <- "screen"
        
        if (!{"arm" %in% colnames(data.clinical)}) {
          arm <- rep("all", nrow(data.clinical))
        } else {
          arm <- data.clinical$arm
        }
        
      } else if (all(data.assess$screen_type == data.assess$screen_type[1L])) {
        # if screen-type designations provided and all are the same, use as
        # name for sensitivity variables
        names(theta$beta) <- data.assess$screen_type[1L]
        names(prior$a_beta) <- data.assess$screen_type[1L]
        names(prior$b_beta) <- data.assess$screen_type[1L]
        screen_type <- data.assess$screen_type
        
        if (!{"arm" %in% colnames(data.clinical)}) {
          arm <- rep("all", nrow(data.clinical))
        } else {
          arm <- data.clinical$arm
        }
        
      } else {
        # if screen-type designations provided and all are not the same, warn
        # the user and reset to a default value
        message("*** only 1 sensitivity type specified in theta; screen_type ignored ***")
        warning("*** only 1 sensitivity type specified in theta; screen_type ignored ***",
                call. = FALSE)
        # reset to a default value
        screen_type <- rep("screen", nrow(data.assess))
        names(theta$beta) <- "screen"
        names(prior$a_beta) <- "screen"
        names(prior$b_beta) <- "screen"
        
        if (!{"arm" %in% colnames(data.clinical)}) {
          arm <- rep("all", nrow(data.clinical))
        } else {
          arm <- data.clinical$arm
        }
      }
    } else {
      # if theta$beta is a vector numeric, multiple arms or screen types. Names
      # can indicate either the study arm or the screen type
      if ({"arm" %in% colnames(data.clinical)} && all(names(theta$beta) %in% data.clinical$arm)) {
        # theta$beta pertains to study arm designation -- assign same screen
        # designation to all screens of each arm
        message("using arm specific sensitivity parameters")

        if ({"screen_type" %in% colnames(data.assess)} &&
            !all(data.assess$screen_type == data.assess$screen_type[1L])) {
          message("*** screen_type designations ignored ***")
          warning('*** screen_type designations ignored ***', call. = FALSE)
        }
        
        screen_type <- rep(NA_character_, nrow(data.assess))
        for (i in 1L:length(theta$beta)) {
          group <- data.clinical$arm == names(theta$beta)[i]
          screen_type[data.assess$id %in% data.clinical$id[group]] <- names(theta$beta)[i]
        }
        if (any(is.na(screen_type))) {
          stop("some screens do not fall into defined arm designations", call. = FALSE)
        }
        arm <- data.clinical$arm
      } else if ({"screen_type" %in% colnames(data.assess)} && 
                 all(data.assess$screen_type %in% names(theta$beta))) {
        # theta pertains to screening type -- assign arm if none given. 
        screen_type <- data.assess$screen_type
        if (!{"arm" %in% colnames(data.clinical)}) {
          arm <- rep("all", nrow(data.clinical))
        } else {
          arm <- data.clinical$arm
        }
      } else {
        stop("could not connect theta$beta to one of data.clinical$arm or data.assess$screen_type; verify input",
             call. = FALSE)
      }
    }
  } else {
    stop("theta$beta must be a numeric vector", call. = FALSE)
  }
  
  list("screen.type" = screen_type, 
       "arm" = arm, 
       "theta" = theta, 
       "prior" = prior)
}



.testGroupBetaInput <- function(screen.type, theta, prior) {
  
  # convert screen type variables to factor and identify unique types
  screen.type <- factor(screen.type)
  all_screen_types <- levels(screen.type)
  
  # ensure that prior and theta parameters are all of the same length
  if (length(prior$a_beta) != length(theta$beta) ||
      length(prior$b_beta) != length(theta$beta)) {
    stop("inappropriate parameter specifications for beta; ",
         "verify `theta_0$beta`, `prior$a_beta`, and `prior$b_beta`", 
         call. = FALSE)
  }
  
  # ensure that all theta and prior vectors have the same length and include
  # all of the types of screens
  if (any(!{names(theta$beta) %in% all_screen_types}) ||
      any(!{names(prior$a_beta) %in% all_screen_types}) ||
      any(!{names(prior$b_beta) %in% all_screen_types}) ||
      any(!{all_screen_types %in% names(theta$beta)}) ||
      any(!{all_screen_types %in% names(prior$a_beta)}) ||
      any(!{all_screen_types %in% names(prior$b_beta)})) {
    stop("names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` ",
         "do not match values of `screen_type`\n\t", 
         "from  `screen_type`: ", paste(all_screen_types, collapse = ", "), "\n\t",
         "from `theta_0$beta`: ", paste(names(theta$beta), collapse = ", "), "\n\t",
         "from `prior$a_beta`: ", paste(names(prior$a_beta), collapse = ", "), "\n\t",
         "from `prior$b_beta`: ", paste(names(prior$b_beta), collapse = ", "), "\n\t",
         call. = FALSE)
  }
  
  # reorder theta and prior vectors to match factor levels
  theta$beta <- theta$beta[all_screen_types]
  prior$a_beta <- prior$a_beta[all_screen_types]
  prior$b_beta <- prior$b_beta[all_screen_types]
  
  list("screen.type" = screen.type,
       "theta" = theta, 
       "prior" = prior)
}

.testGroupRatePInput <- function(grp.rateP, n, theta, prior, epsilon = NULL) {
  
  # to handle calls from functions that do not perform the mcmc step
  if (is.null(epsilon)) epsilon <- theta$rate_P
  
  # if no grp.rateP specified in data.clinical, set to default
  if (is.null(grp.rateP)) {
    message("no grouping variable for pre-clinical Weibull distribution specified\n\t",
            "all participants assumed to follow the same pre-clinical distribution")
    grp.rateP <- rep("all", n)
  }
  
  if (!is.factor(grp.rateP)) grp.rateP <- factor(grp.rateP)
  
  groups_rateP <- levels(grp.rateP)
  
  # if theta variable has only 1 value, but more than 1 rateP group in the
  # data, warn user and reset rateP grouping variable
  if (length(theta$rate_P) == 1L && length(groups_rateP) != 1L) {
    message("`grp.rate_P` variable of `data.clinical` ignored; ",
            "rate_P is not defined as group-specific in theta")
    warning("`grp.rate_P` variable of `data.clinical` ignored; ",
            "rate_P is not defined as group-specific in theta\n", call. = FALSE)
    grp.rateP <- factor(rep("all", n))
    groups_rateP <- levels(grp.rateP)
  }

  # ensure that all rateP related parameters are of the same length
  if (length(prior$rate_P) != length(theta$rate_P) ||
      length(prior$shape_P) != length(theta$rate_P) ||
      length(epsilon) != length(theta$rate_P)) {
    stop("inappropriate parameter specifications for pre-clinical Weibull; ",
         "verify `theta_0$rate_P`, `prior$rate_P`, prior$shape_P, and `epsilon_rate_P` ",
         call. = FALSE)
  }
  
  if (length(theta$rate_P) > 1L) {
    # if more than 1 rate_p provided, verify that all names match
    if (any(!{names(theta$rate_P) %in% groups_rateP}) ||
        any(!{names(prior$rate_P) %in% groups_rateP}) ||
        any(!{names(prior$shape_P) %in% groups_rateP}) ||
        any(!{names(epsilon) %in% groups_rateP}) ||
        any(!{groups_rateP %in% names(theta$rate_P)}) ||
        any(!{groups_rateP %in% names(prior$rate_P)}) ||
        any(!{groups_rateP %in% names(prior$shape_P)}) ||
        any(!{groups_rateP %in% names(epsilon)})) {
      stop("names of `theta_0$rate_P`, `prior$rate_P`, `prior$shape_P`, and/or `epsilon_rate_P` ",
           "do not match values of `grp.rateP`\n\t", 
           "from `grp.rateP`: ", paste(groups_rateP, collapse = ", "), "\n\t",
           "from `theta_0$rate_P`: ", paste(names(theta$rate_P), collapse = ", "), "\n\t",
           "from `prior$rate_P`: ", paste(names(prior$rate_P), collapse = ", "), "\n\t",
           "from `prior$shape_P`: ", paste(names(prior$shape_P), collapse = ", "), "\n\t",
           "from `epsilon_rate_P`: ", paste(names(epsilon), collapse = ", "), "\n\t", 
           call. = FALSE)
    }
    theta$rate_P <- theta$rate_P[groups_rateP]
    prior$rate_P <- prior$rate_P[groups_rateP]
    prior$shape_P <- prior$shape_P[groups_rateP]
    epsilon <- epsilon[groups_rateP]
  } else {
    names(theta$rate_P) <- groups_rateP
    names(prior$rate_P) <- groups_rateP
    names(prior$shape_P) <- groups_rateP
    names(epsilon) <- groups_rateP
  }
  
  list("grp.rateP" = grp.rateP, 
       "theta" = theta, 
       "prior" = prior, 
       "epsilon" = epsilon)
}

.testAdaptiveInput <- function(adaptive) {
  
  if (!is.vector(adaptive, "list")) {
    stop("`adaptive` must be null or a list", call. = FALSE)
  }
  
  if (!all(names(adaptive) %in% c("warmup", "delta", "gamma", "kappa", "m0")) ||
      !all(c("warmup", "delta", "gamma", "kappa", "m0") %in% names(adaptive))) {
    stop("`adaptive` does not contain required information", call. = FALSE)
  }
  
  if (any(!{lapply(adaptive, is.vector, mode = "numeric") |> unlist()})) {
    stop("all elements of `adaptive` must be numeric", call. = FALSE)
  }
  
  if (any(lengths(adaptive) > 1L)) {
    stop("all elements of `adaptive` must be scalar", call. = FALSE)
  }
  
  if (!isTRUE(all.equal(adaptive$warmup, as.integer(adaptive$warmup))) ||
      adaptive$warmup < 0) {
    stop("`adaptive$warmup` must be a positive integer", call. = FALSE)
  }
  
  if (adaptive$delta < 0 || adaptive$delta > 1.0) {
    stop("`adaptive$delta` must in (0, 1)", call. = FALSE)
  }
  
  if (adaptive$gamma < 0 || adaptive$gamma > 1.0) {
    stop("`adaptive$gamma` must in (0, 1]", call. = FALSE)
  }
  
  if (adaptive$kappa < 0 || adaptive$kappa > 1.0) {
    stop("`adaptive$kappa` must in (0, 1]", call. = FALSE)
  }
  
  if (!isTRUE(all.equal(adaptive$m0, as.integer(adaptive$m0))) ||
      adaptive$m0 < 0) {
    stop("`adaptive$m0` must be a positive integer", call. = FALSE)
  }
  
}


.testDataAssess <- function(data.assess) {
  
  data.assess$id <- data.assess$id |> as.character()

  if (!is.numeric(data.assess$disease_detected)) {
    stop("data.assess$disease_detected must be binary 0/1", call. = FALSE)
  }
  
  if (!is.integer(data.assess$disease_detected)) {
    data.assess$disease_detected <- as.integer(round(data.assess$disease_detected))
  }

  if (!all(data.assess$disease_detected %in% c(0L, 1L))) {
    stop("`data.assess$disease_detected` must be binary 0/1", call. = FALSE)
  }
  
  # ages cannot be negative
  if (any(data.assess$age_assess < 0.0)) {
    stop("`data.assess$age_assess` cannot be negative", call. = FALSE)
  }
  
  # order the history and clinical datasets
  data.assess <- data.assess[order(data.assess$id, data.assess$age_assess), ]
  rownames(data.assess) <- NULL
  data.assess
}

.testDataClinical <- function(data.clinical) {
  
  data.clinical$id <- data.clinical$id |> as.character()
  
  # ensure only 1 record per id is present in clinical data
  if (any(duplicated(data.clinical$id))) {
    stop("at least 1 participant id has multiple records in `data.clinical`",
         call. = FALSE)
  }
  
  # ages cannot be negative
  if (any(data.clinical$age_endpoint < 0.0)) {
    stop("`data.clinical$age_endpoint` cannot be negative", call. = FALSE)
  }
  
  if (any(data.clinical$age_entry < 0.0)) {
    stop("`data.clinical$age_entry` cannot be negative", call. = FALSE)
  }
  
  # ensure all endpoint types are as expected
  if (!all(data.clinical$endpoint_type %in% c("clinical", "preclinical", "censored"))) {
    stop("unrecognized endpoint type in `data.clinical$endpoint_type`", call. = FALSE)
  }
  
  # order the history and clinical datasets
  data.clinical <- data.clinical[order(data.clinical$id), ]
  rownames(data.clinical) <- NULL
  data.clinical
}

.testThetaInput <- function(theta) {
  
  # most elements must be scalar numeric
  if (any(lengths(theta[c("rate_H", "shape_H", "shape_P", "psi")]) > 1L)) {
    stop("`theta_0` elements rate_H, shape_H, shape_P, and psi must be scalars", 
         call. = FALSE)
  }
  
  tst <- lapply(theta, function(x) { all(x >= 0.0) }) |> unlist()
  
  if (any(!tst)) {
    stop("initial values for distribution parameters must be non-negative; ", 
         "verify `theta_0`", call. = FALSE)
  }
  
  if (any(theta$beta > 1.0)) {
    stop("`theta_0$beta` must be in [0, 1]", call. = FALSE)
  }
  if (theta$psi > 1.0) stop("`theta_0$psi` must be in [0, 1]", call. = FALSE)
  
  theta
}

.testPriorInput <- function(prior) {
  # most elements must be scalar numeric
  if (any(lengths(prior[c("rate_H", "shape_H", "a_psi", "b_psi")]) > 1L)) {
    stop("`prior` elements rate_H, shape_H, a_psi, and b_psi must be scalars", 
         call. = FALSE)
  }
  
  if (length(prior$shape_P) != length(prior$rate_P)) {
    stop("lengths of prior$shape_P and prior$rate_P must be the same", call. = FALSE)
  }
  
  if (length(prior$a_beta) != length(prior$b_beta)) {
    stop("lengths of prior$a_beta and prior$b_beta must be the same", call. = FALSE)
  }

  tst <- lapply(prior, function(x) { all(x >= 0.0) }) |> unlist()
  
  if (any(!tst)) {
    stop("all distribution parameters of the priors must be non-negative; ", 
         "verify `prior`", call. = FALSE)
  }
  
  prior
}

.getPreviousResults <- function(data.objects, baclava.object, indolent) {
  
  cnt <- 1L
  
  age_at_tau_hats <- list()
  indolents <- list()
  for (obj in data.objects) {
    idx_obj <- match(obj$ids, colnames(baclava.object$tau_hp))
    if (any(is.na(idx_obj))) stop("data ids do not match previous result", call. = FALSE)
    
    last_row <- nrow(baclava.object$tau_hp)
    while (TRUE) {
      age_at_tau_hats[[cnt]] <- baclava.object$tau_hp[last_row, idx_obj]
      if (!all(age_at_tau_hats[[cnt]] < 1e-8)) break
      last_row <- last_row - 1L
      if (last_row <= 0L) stop("all prior tau_hp values are zero", call. = FALSE)
    }
    
    if (indolent) indolents[[cnt]] <- baclava.object$indolent[last_row, idx_obj]
    
    cnt <- cnt + 1L
  }
  
  list("age.at.tau.hp.hats" = age_at_tau_hats, "indolents" = indolents)      
  
}