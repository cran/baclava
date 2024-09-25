# add subscript to rate_P and beta related names when appropriate
.getNames <- function(theta) {
  nms <- list(
    "beta" = if (ncol(theta$beta) > 1L) {
               paste0("beta[", colnames(theta$beta), "]") 
             } else { "beta" },
    "psi" = "psi", 
    "lambda[H]" = "lambda[H]", 
    "mu[H]" = "mu[H]",
    "lambda[P]" = if (ncol(theta$"lambda[P]") > 1L) {
                    paste0("lambda[P][", colnames(theta$"lambda[P]"), "]") 
                  } else { "lambda[P]" }, 
    "mu[P]" = if (ncol(theta$"mu[P]") > 1L) {
                paste0("mu[P][", colnames(theta$"mu[P]"), "]") 
              } else { "mu[P]" }
    )
  nms[names(theta)] |> unlist()
}

# convert coding names to strings provided for expression()
.makePrettyNames <- function(theta) {
  pretty_names <- c("beta" = "beta", "psi" = "psi", "rate_H" = "lambda[H]",
                    "rate_P" = "lambda[P]", "mean_H" = "mu[H]",
                    "mean_P" = "mu[P]")
  
  # reset names of theta list to be "pretty"
  idx <- match(names(theta), names(pretty_names))
  names(theta) <- pretty_names[idx]
  theta
}

#'
#' @describeIn fit_baclava Summary statistics of posterior distribution parameters
#' @param object An object of class \code{baclava}.
#' @param ... Ignored.
#' @importFrom coda effectiveSize
#' @importFrom stats quantile
#' @importFrom tibble as_tibble
#' @export
summary.baclava <- function(object, ...) {
  
  object$theta$mean_H <- object$theta$rate_H^{-1.0 / object$theta$shape_H} *
    gamma(1.0 + 1.0 / object$theta$shape_H)
  
  object$theta$mean_P <- object$theta$rate_P^{-1.0 / object$theta$shape_P} *
    gamma(1.0 + 1.0 / object$theta$shape_P)
  
  # we don't care about shape parameters here
  object$theta$shape_H <- NULL
  object$theta$shape_P <- NULL

  # convert names from coding choices to values useful in expression()
  object$theta <- .makePrettyNames(object$theta)

  # if disease always progressive, get rid of psi
  if (!object$setup$indolent) x$theta$psi <- NULL
  
  .summary <- function(x) {
    c("mean" = mean(x),
      "q_low" = stats::quantile(x, probs = 0.025),
      "q_high" = stats::quantile(x, probs = 0.975),
      "ESS" = coda::effectiveSize(x))
  }
  
  res <- lapply(object$theta, 
                function(x) {
                  if (is.matrix(x)) {
                    apply(x, 2L, .summary)
                  } else {
                    .summary(x)
                  }
                }) |> do.call(what = "cbind")

  res <- res |> 
    matrix(ncol = 4L, byrow = TRUE,
           dimnames = list(.getNames(object$theta), c("mean", "q_low", "q_high", "ESS")))

  tibble::as_tibble(res, rownames = "param")

}

#'
#' @describeIn fit_baclava Print summary statistics of posterior distribution parameters
#' @param x An object of class \code{baclava}.
#' @export
print.baclava <- function(x, ...) {
  print(summary.baclava(object = x))
  invisible(x)
}

#' Plot Posterior Distribution Parameters
#'
#' Convenience function to facilitate exploration of posterior distributions through
#'   trace plots, autocorrelations, and densities, as well as plotting the 
#'   estimated hazard for transitioning to the preclinical compartment.
#'
#' @param x An object of class \code{baclava}.
#' @param y Ignored
#' @param type A character object. One of \{"density", "trace", "acf", "hazard"\}. The
#'   type of plot to generate
#' @param burnin An integer object. Optional. The number of burn-in samples.
#'   Used only for \code{type = "trace"}. One trace plot is generated for
#'   the burnin iterations; a second for the post-burnin iterations. Note, this
#'   refers to the number of kept (thinned) samples.
#' @param max_age A numeric object. For \code{type = "hazard"}, the maximum
#'   age at which to evaluate the hazard.
#' @param trace_var A character object. The parameter for which trace plots
#'   are to be generated. Must be one of \{"psi", "rate_H", "rate_P", "beta"\}
#' @param ... Ignored
#' 
#' @returns A gg object
#' 
#' @examples
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
#' plot(example)
#' plot(example, type = "trace", trace_var = "psi", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "rate_H", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "rate_P", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "beta", burnin = 0L) 
#' plot(example, type = "acf")
#' plot(example, type = "hazard", max_age = 70)
#' 
#' @import ggplot2
#' @importFrom stats acf quantile
#' @export
plot.baclava <- function(x, y, ..., 
                         type = c("density", "trace", "acf", "hazard"), 
                         burnin = 0L,
                         max_age = 90L,
                         trace_var = c("psi", "rate_H", "rate_P", "beta")) {

  type <- match.arg(type)
  
  x$theta$mean_H <- x$theta$rate_H^{-1.0 / x$theta$shape_H} *
    gamma(1.0 + 1.0 / x$theta$shape_H)
  
  x$theta$mean_P <- matrix(x$theta$rate_P^{-1.0 / x$theta$shape_P} *
                             gamma(1.0 + 1.0 / x$theta$shape_P), 
                           ncol = ncol(x$theta$rate_P),
                           dimnames = list(NULL, colnames(x$theta$rate_P)))
  
  shape_H <- x$theta$shape_H
  shape_P <- x$theta$shape_P
  
  pretty_names <- c("beta" = "beta", "psi" = "psi", 
                    "rate_H" = "lambda[H]", "mean_H" = "mu[H]", 
                    "rate_P" = "lambda[P]", "mean_P" = "mu[P]")
  x$theta <- x$theta[names(pretty_names)]
  if (!x$setup$indolent) x$theta$psi <- NULL

  if (type == "density") {

    # extract all elements as a data.frame
    res <- lapply(1L:length(x$theta), 
                  function(i, theta) {
                    if (is.matrix(theta[[i]])) {
                      data.frame("value" = c(theta[[i]]),
                                 "parm" = rep(paste0(pretty_names[i], "['", 
                                                     colnames(theta[[i]]), "']"),
                                              each = nrow(theta[[i]])))
                    } else {
                      data.frame("value" = theta[[i]], 
                                 "parm" = unname(pretty_names[i]))
                    }
                  },
                  theta = x$theta)

    df <- do.call(rbind, res)
    df$parm <- factor(df$parm)

    df_summary_names <- sapply(1L:length(x$theta), 
                               function(i, theta) {
                                 if (is.matrix(theta[[i]])) {
                                   paste0(pretty_names[i], "['", 
                                          colnames(theta[[i]]), "']")
                                 } else {
                                   unname(pretty_names[i])
                                 }
                               },
                               theta = x$theta) |> unlist()
    
    # get summary statistics for each parameter
    df_summary <- lapply(x$theta, 
                         function(x) {
                           if (is.matrix(x)) { colMeans(x) } else { mean(x) }
                         }) |> unlist()

    df_summary <- data.frame("mean" = df_summary, "parm" = factor(df_summary_names))

    df_summary$Q1 <- lapply(x$theta, 
                            function(tmp) {
                              if (is.matrix(tmp)) {
                                apply(tmp, 2L, stats::quantile, probs = 0.025)
                              } else {
                                stats::quantile(tmp, probs = 0.025)
                              }
                            }) |> unlist()

    df_summary$Q2 <- lapply(x$theta, 
                            function(tmp) {
                              if (is.matrix(tmp)) {
                                apply(tmp, 2L, stats::quantile, probs = 0.975)
                              } else {
                                stats::quantile(tmp, probs = 0.975)
                              }
                            }) |> unlist()
    
    ggplot(data = df) +
      geom_density(aes(.data$value)) +
      facet_wrap(.~parm, labeller = label_parsed, scales = "free", ncol = 2) +
      labs(x = "") +
      geom_vline(data = df_summary, aes(xintercept = mean, group = .data$parm), colour = "maroon") +
      geom_vline(data = df_summary, aes(xintercept = .data$Q1, group = .data$parm), colour = "maroon") +
      geom_vline(data = df_summary, aes(xintercept = .data$Q2, group = .data$parm), colour = "maroon")
    
  } else if (type == "trace") {
    
    trace_var <- match.arg(trace_var)
    idx <- match(trace_var, names(x$theta))

    if (is.na(idx)) {
      stop("requested parameter is not present in `$theta`", call. = FALSE)
    }
    
    if (burnin < 0 || 
        burnin > ifelse(is.matrix(x$theta[[idx]]), 
                        nrow(x$theta[[idx]]), 
                        length(x$theta[[idx]]))) {
      stop("burnin is not appropriately set", call. = FALSE)
    }
    
    res <- lapply(idx, 
                  function(i, theta) {
                    if (is.matrix(theta[[i]])) {
                      data.frame("sample" = seq_len(nrow(theta[[i]])),
                                 "value" = c(theta[[i]]),
                                 "parm" = rep(paste0(pretty_names[i], "['", 
                                                     colnames(theta[[i]]), "']"),
                                              each = nrow(theta[[i]])))
                    } else {
                      data.frame("sample" = seq_along(theta[[i]]),
                                 "value" = theta[[i]], 
                                 "parm" = unname(pretty_names[i]))
                    }
                  },
                  theta = x$theta)
    df <- do.call(rbind, res)
    df$parm <- factor(df$parm)
    
    df$group <- factor(c("Recurrent", "Transient")[as.integer(df$sample <= burnin) + 1L],
                       levels = c("Transient", "Recurrent"))

    ggplot(df) + 
      geom_line(aes(x = .data$sample, y = .data$value)) +
      labs(x = "Sample", y = "") +
      scale_x_continuous(trans = ifelse(nrow(df) > 1000, "log10", "identity")) +
      facet_wrap(~group + parm, labeller = label_parsed, scales = "free_x")

  } else if (type == "acf") {
    
    my_acf <- function(x) {
      tmp <- stats::acf(x = x, plot = FALSE)
      data.frame("acf" = as.vector(tmp$acf), "lag" = as.vector(tmp$lag))
    }
    
    df <- lapply(seq_along(x$theta), 
                 function(i, theta) {
                   if (is.matrix(theta[[i]])) {
                     res <- apply(theta[[i]], 2, my_acf, simplify = FALSE)
                     res2 <- do.call(rbind, res)
                     res2$parm <- rep(paste0(pretty_names[i], "['", 
                                             colnames(theta[[i]]), "']"),
                                      each = nrow(res[[1L]]))
                     res2
                   } else {
                     res2 <- my_acf(theta[[i]])
                     res2$parm <- pretty_names[i]
                     res2
                   }
                 },
                 theta = x$theta) 
    df <- do.call(rbind, df)
    df$parm <- factor(df$parm)
    
    ggplot() + geom_col(data = df, aes(.data$lag, .data$acf)) + 
      scale_y_continuous(breaks = seq(0.0, 1.0, 0.5)) +
      facet_wrap(~parm, labeller = label_parsed, ncol = 2L) +
      labs(x = "Lags", y= "Auto-correlation\n function")
    
  } else if (type == "hazard") {
    
    if (max_age <= x$setup$t0) {
      stop("max_age must be greater than the risk onset age ", x$t0, call. = FALSE)
    }
    rate_H_quantile <- stats::quantile(x$theta$rate_H, c(0.025, 0.5, 0.975))
    df <- data.frame("ages" = seq(x$setup$t0, max_age, by = 0.1))
    tmp <- {df$ages - x$setup$t0}^{shape_H - 1L}
    df$h_low <- shape_H * rate_H_quantile[1L] * tmp
    df$h_med <- shape_H * rate_H_quantile[2L] * tmp
    df$h_upp <- shape_H * rate_H_quantile[3L] * tmp
    
    ggplot(data = df, aes(x = .data$ages)) +
      geom_line(aes(y = .data$h_med)) +
      geom_line(aes(y = .data$h_low), linetype = "dashed") +
      geom_line(aes(y = .data$h_upp), linetype = "dashed") +
      labs(x="Age", y = "Hazard rate for\npre-clinical cancer")
  }
}
