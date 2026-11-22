#' @title bdpnormal Object Summary
#' @description \code{summary} method for class \code{bdpnormal}.
#' @param object object of class \code{bdpnormal}. The result of a call to the
#'   \code{\link{bdpnormal}} function.
#'
#' @details Displays input data, the current-historical comparison, and
#'   historical weights when both sources are available. Under \code{method =
#'   "mc"}, the comparison and weight summaries are medians across draws.
#'   The outcome estimates are posterior medians of each arm's population mean.
#'   Intervals are the 2.5th and 97.5th posterior percentiles. For a two-arm
#'   analysis, the interval describes the treatment-minus-control difference,
#'   calculated from paired draws. It is not a confidence interval or a test
#'   of equality. Printed estimates retain four decimal places.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats sd density is.empty.model median model.offset
#'   model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("summary", signature(object = "bdpnormal"), function(object) {
  estimates <- bdp_arm_summary(object, "posterior_mu")
  arm2 <- object$args1$arm2
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  mu_t <- object$args1$mu_t
  sigma_t <- object$args1$sigma_t
  N_t <- object$args1$N_t
  mu0_t <- object$args1$mu0_t
  sigma0_t <- object$args1$sigma0_t
  N0_t <- object$args1$N0_t
  mu_c <- object$args1$mu_c
  sigma_c <- object$args1$sigma_c
  N_c <- object$args1$N_c
  mu0_c <- object$args1$mu0_c
  sigma0_c <- object$args1$sigma0_c
  N0_c <- object$args1$N0_c
  method <- object$args1$method
  
  if (!arm2) {
    # Format treatment mean output
    mean_CI_t <- round(c(estimates$conf.low[1], estimates$conf.high[1]), 4)
    mean_est_t <- round(estimates$estimate[1], 4)
    
    cat("\n")
    cat("    One-armed bdp normal\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: mu_t = ", mu_t, ", sigma_t = ", sigma_t, ", N_t = ", N_t))
    cat("\n")
    if (!is.null(N0_t)) {
      cat(paste0(
        "  Historical treatment: mu0_t = ", mu0_t, ", sigma0_t = ",
        sigma0_t, ", N0_t = ", N0_t
      ))
      cat("\n")
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
    }
    cat("95 percent equal-tailed credible interval: \n")
    cat(paste0(" ", mean_CI_t))
    cat("\n")
    cat("Posterior median:\n")
    cat("population mean of treatment group\n")
    cat(paste0(" ", mean_est_t))
    cat("\n")
  } else {
    mean_est_t <- format(round(estimates$estimate[1], 4), nsmall = 4)
    mean_est_c <- format(round(estimates$estimate[2], 4), nsmall = 4)
    comp_CI <- round(c(estimates$conf.low[3], estimates$conf.high[3]), 4)
    cat("\n")
    cat("    Two-armed bdp normal\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: mu_t = ", mu_t, ", sigma_t = ", sigma_t, ", N_t = ", N_t))
    cat("\n")
    if (!is.null(N_c)) {
      cat(paste0("  Current control: mu_c = ", mu_c, ", sigma_c = ", sigma_c, ", N_c = ", N_c))
      cat("\n")
    }
    if (!is.null(N0_t)) {
      cat(paste0("  Historical treatment: mu0_t = ", mu0_t, ", sigma0_t = ", sigma0_t, ", N0_t = ", N0_t))
      cat("\n")
    }
    if (!is.null(N0_c)) {
      cat(paste0("  Historical control: mu0_c = ", mu0_c, ", sigma0_c = ", sigma0_c, ", N0_c = ", N0_c))
      cat("\n")
    }
    
    if (!is.null(N0_t)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(N0_c) & !is.null(N_c)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - control (current vs. historical data): ",
        round(median(posterior_control$p_hat), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(N0_t)) {
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(N0_c) & !is.null(N_c)) {
      cat(paste0(
        "Discount function value (alpha) - control: ",
        round(median(posterior_control$alpha_discount), 4)
      ))
      cat("\n")
    }
    cat("Treatment effect: treatment - control")
    cat("\n")
    
    cat("95 percent equal-tailed credible interval: \n")
    cat(paste0(" ", comp_CI))
    cat("\n")
    cat("Posterior medians:\n")
    cat("treatment group  control group\n")
    cat(paste0("          ", mean_est_t, "          ", mean_est_c))
    cat("\n")
  }
})


#' @title bdpbinomial Object Summary
#' @description \code{summary} method for class \code{bdpbinomial}.
#' @param object object of class \code{bdpbinomial}. The result of a call to the
#'   \code{\link{bdpbinomial}} function.
#'
#' @details Displays input data, the current-historical comparison, and
#'   historical weights when both sources are available. Under \code{method =
#'   "mc"}, the comparison and weight summaries are medians across draws.
#'   The outcome estimates are posterior medians of each arm's event probability.
#'   Intervals are the 2.5th and 97.5th posterior percentiles. For a two-arm
#'   analysis, the interval describes the treatment-minus-control difference,
#'   calculated from paired draws. It is not a confidence interval or a test
#'   of equality. Printed estimates retain four decimal places.
#'
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull quantile rbeta rgamma rnorm var vcov
#' @export
setMethod("summary", signature(object = "bdpbinomial"), function(object) {
  estimates <- bdp_arm_summary(object, "posterior")
  arm2 <- object$args1$arm2
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  y_t <- object$args1$y_t
  N_t <- object$args1$N_t
  y_c <- object$args1$y_c
  N_c <- object$args1$N_c
  y0_t <- object$args1$y0_t
  N0_t <- object$args1$N0_t
  y0_c <- object$args1$y0_c
  N0_c <- object$args1$N0_c
  method <- object$args1$method
  
  if (!arm2) {
    ### Format treatment mean output
    mean_CI_t <- round(c(estimates$conf.low[1], estimates$conf.high[1]), 4)
    mean_est_t <- round(estimates$estimate[1], 4)
    
    cat("\n")
    cat("    One-armed bdp binomial\n\n")
    cat(paste0("Current treatment data: ", y_t, " and ", N_t, "\n"))
    if (!is.null(N0_t)) {
      cat(paste0("Historical treatment data: ", y0_t, " and ", N0_t, "\n"))
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
    }
    cat("95 percent equal-tailed credible interval: \n")
    cat(paste0(" ", mean_CI_t))
    cat("\n")
    cat("Posterior medians:\n")
    cat(paste0(" ", mean_est_t))
    cat("\n")
  } else {
    
    ### Format control mean output
    mean_est_t <- format(round(estimates$estimate[1], 4), nsmall = 4)
    mean_est_c <- format(round(estimates$estimate[2], 4), nsmall = 4)
    
    ### Format comparison output
    comp_CI <- round(c(estimates$conf.low[3], estimates$conf.high[3]), 4)
    
    cat("\n")
    cat("    Two-armed bdp binomial\n\n")
    cat(paste0("Current treatment data: ", y_t, " and ", N_t, "\n"))
    if (!is.null(N_c)) {
      cat(paste0("Current control data: ", y_c, " and ", N_c, "\n"))
    }
    if (!is.null(N0_t)) {
      cat(paste0("Historical treatment data: ", y0_t, " and ", N0_t, "\n"))
    }
    if (!is.null(N0_c)) {
      cat(paste0("Historical control data: ", y0_c, " and ", N0_c, "\n"))
    }
    
    if (!is.null(N_t) & !is.null(N0_t)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(N_c) & !is.null(N0_c)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - control (current vs. historical data): ",
        round(median(posterior_control$p_hat), 4)
      ))
      cat("\n")
      cat(paste0(
        "Discount function value (alpha) - control: ",
        round(median(posterior_control$alpha_discount), 4)
      ))
      cat("\n")
    }
    
    cat("Treatment effect: treatment - control")
    cat("\n")
    cat("95 percent equal-tailed credible interval:\n")
    cat(paste0(" "), comp_CI)
    cat("\n")
    cat("Posterior medians:\n")
    cat("treatment  control\n")
    cat(paste0("  ", mean_est_t, "  ", mean_est_c))
  }
})


#' @title bdpsurvival Object Summary
#' @description \code{summary} method for class \code{bdpsurvival}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull quantile rbeta rgamma rnorm var vcov
#' @param object an object of class \code{bdpsurvival}, a result of a call to
#'   \code{\link{bdpsurvival}}.
#' @details Displays a summary of the \code{bdpsurvival} fit. The output is
#'   different, conditional on a one- or two-arm survival analysis.
#'
#'   In the case of a one-arm analysis, the stochastic comparison between
#'   current and historical data and the resulting historical data weight
#'   (alpha) are displayed, followed by a survival table containing augmented
#'   posterior estimates of the survival probability at each time point for the
#'   current data. Times are observed event or censoring times. Risk sets and
#'   event counts describe current participants, not the split hazard intervals.
#'   The survival estimates are posterior medians, with posterior SDs and 95%
#'   equal-tailed credible intervals; they are not Kaplan--Meier estimates.
#'   Saved fits without participant-level records must be refitted before
#'   calling \code{summary()}.
#'
#'   When a control arm is present, the summary describes the approximate
#'   common log hazard ratio from the piecewise hazard draws. The comparison
#'   between current and historical data and the resulting historical data
#'   weight (alpha) are displayed, when historical data is present for the
#'   respective arm. The displayed \code{coef} is the posterior mean log hazard
#'   ratio, treatment relative to control; \code{exp(coef)} exponentiates this
#'   mean, and \code{se(coef)} is its posterior SD. The 95 percent equal-tailed
#'   credible limits are on the log hazard ratio scale. Sample sizes and event
#'   counts describe the current arms, or historical controls when current
#'   controls are absent.
#'
#' @export
setMethod("summary", signature(object = "bdpsurvival"), function(object) {
  posterior_treatment <- object$posterior_treatment
  posterior_control <- object$posterior_control
  surv_time <- object$args1$surv_time
  
  args1 <- object$args1
  data <- args1$participants
  if (is.null(data)) {
    stop("This fit lacks participant-level records. Refit with the current version of bayesDP before calling summary().", call. = FALSE)
  }
  breaks <- args1$breaks
  arm2 <- args1$arm2
  method <- args1$method
  
  historical <- NULL
  treatment <- NULL
  
  ##############################################################################
  ### Survival table
  ### - Only print if !arm2
  ##############################################################################
  if (!arm2) {
    ### Organize current treatment posterior data
    data_t <- subset(data, historical == 0 & treatment == 1)
    s_t <- survival::survfit(Surv(time, status) ~ 1, data = data_t,
                            timefix = FALSE)
    
    estimates <- bdp_survival_summary(object, s_t$time)
    s_t$surv <- estimates$estimate
    s_t$std.err <- estimates$std.error
    s_t$upper <- estimates$conf.high
    s_t$lower <- estimates$conf.low
    
    m_t <- round(cbind(
      s_t$time, s_t$n.risk, s_t$n.event, s_t$surv, s_t$std.err,
      s_t$lower, s_t$upper
    ), 4)
    cnames <- c(
      "time", "n.risk", "n.event", "median", "SD",
      "lower 95% CrI", "upper 95% CrI"
    )
    dimnames(m_t) <- list(rep("", nrow(m_t)), cnames)
    
    cat("\n")
    cat("    One-armed bdp survival\n\n")
    cat("Posterior median survival, SD, and 95% equal-tailed credible interval\n")
    if (is.null(args1$S0_t)) {
      cat("  Current treatment summary:")
      cat("\n")
      print(m_t)
    } else {
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
      
      cat("\n")
      cat("Current treatment - augmented posterior summary:")
      cat("\n")
      print(m_t)
      cat("\n")
    }
  }
  
  ##############################################################################
  ### Approximate common log hazard ratio for treatment vs control
  ### - Only print if arm2
  ##############################################################################
  if (arm2) {
    ### Compute treatment effect of treatment vs control via the same
    ### inverse-variance-weighted piecewise log-hazard-ratio estimator used in
    ### bdpsurvival(); see the survival vignette, 'Approximate two-arm effect'.
    estimates <- bdp_survival_effect_summary(object)
    coef <- estimates$estimate
    se_coef <- estimates$std.error
    CI <- c(estimates$conf.low, estimates$conf.high)
    
    summ_table <- matrix(round(c(coef, exp(coef), se_coef, CI[1], CI[2]), 4), nrow = 1)
    cnames <- c("coef", "exp(coef)", "se(coef)", "lower 95% CI", "upper 95% CI")
    dimnames(summ_table) <- list("treatment", cnames)
    
    ### Count sample size and number of events
    data_t <- subset(data, historical == 0 & treatment == 1)
    n_t <- nrow(data_t)
    e_t <- sum(data_t$status)
    
    if (!is.null(args1$S_c)) {
      data_c <- subset(data, historical == 0 & treatment == 0)
      n_c <- nrow(data_c)
      e_c <- sum(data_c$status)
    } else {
      data_c <- subset(data, historical == 1 & treatment == 0)
      n_c <- nrow(data_c)
      e_c <- sum(data_c$status)
    }
    
    cat("\n")
    cat("    Two-armed bdp survival\n\n")
    cat("data:\n")
    cat(paste0("  Current treatment: n = ", n_t, ", number of events = ", e_t))
    cat("\n")
    if (!is.null(args1$S_c)) {
      cat(paste0("  Current control: n = ", n_c, ", number of events = ", e_c))
      cat("\n")
    } else {
      cat(paste0("  Historical control: n = ", n_c, ", number of events = ", e_c))
      cat("\n")
    }
    
    if (!is.null(args1$S0_t)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - treatment (current vs. historical data): ",
        round(median(posterior_treatment$p_hat), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(args1$S0_c) & !is.null(args1$S_c)) {
      cat(paste0(
        "Stochastic comparison (p_hat) - control (current vs. historical data): ",
        round(median(posterior_control$p_hat), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(args1$S0_t)) {
      cat(paste0(
        "Discount function value (alpha) - treatment: ",
        round(median(posterior_treatment$alpha_discount), 4)
      ))
      cat("\n")
    }
    
    if (!is.null(args1$S0_c) & !is.null(args1$S_c)) {
      cat(paste0(
        "Discount function value (alpha) - control: ",
        round(median(posterior_control$alpha_discount), 4)
      ))
      cat("\n")
    }
    cat("\n")
    print(summ_table)
  }
})


#' @title bdplm Object Summary
#' @description \code{summary} method for class \code{bdplm}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull quantile rbeta rgamma rnorm var vcov
#' @param object an object of class \code{bdplm}, a result of a call to
#'   \code{\link{bdplm}}.
#' @details Displays a summary of the \code{bdplm} fit. Displayed output is
#'   similar to \code{\link{summary.lm}}. The function call, residuals, and
#'   coefficient table are displayed. The discount function weight estimates are
#'   displayed as well. If \code{method = "mc"}, the median estimate of alpha is
#'   displayed.
#'   Residuals use current outcomes minus predictions from posterior mean
#'   coefficients, with covariates and intercept on the original scale.
#'   The method is also used by \code{bdplogit()}; for those fits, residuals
#'   are observed binary outcomes minus fitted probabilities obtained by
#'   applying the inverse logit to the posterior mean linear predictor.
#'   Saved fits without the covariate centering means must be refitted before
#'   calling \code{summary()}.
#'
#' @export
setMethod("summary", signature(object = "bdplm"), function(object) {
  
  # Format residuals
  covariate_means <- object$args1$covariate_means
  if (is.null(covariate_means)) {
    stop("This fit lacks covariate centering means. Refit with the current version of bayesDP before calling summary().", call. = FALSE)
  }
  data <- object$args1$data
  data <- data[data$historical == 0, ]
  for (nm in names(covariate_means)) {
    data[[nm]] <- data[[nm]] + covariate_means[[nm]]
  }
  coefficient_names <- setdiff(names(object$estimates$coefficients), "sigma")
  X <- as.matrix(data[, coefficient_names, drop = FALSE])
  y <- data$Y
  beta <- unlist(object$estimates$coefficients[coefficient_names])
  fitted <- drop(X %*% beta)
  if (identical(object$args1$family, "binomial")) {
    fitted <- stats::plogis(fitted)
  }
  resid <- y - fitted
  resid_sum <- data.frame(
    Min = min(resid),
    Q1 = quantile(resid, 0.25),
    Median = median(resid),
    Q3 = quantile(resid, 0.75),
    Max = max(resid)
  )
  resid_sum <- round(resid_sum, 3)
  dimnames(resid_sum) <- list(
    rep("", nrow(resid_sum)),
    c("Min", "1Q", "Median", "3Q", "Max")
  )
  
  # Format coefficient table
  estimates <- bdp_regression_summary(object, conf.int = FALSE)
  coefs <- cbind(estimates$estimate, estimates$std.error)
  dimnames(coefs) <- list(estimates$term, c("Estimate", "Std. Error"))
  coefs <- round(coefs, 4)
  
  # Format alpha
  alpha_mat <- apply(object$alpha_discount, 2, median)
  alpha_mat <- matrix(alpha_mat, nrow = 1)
  dimnames(alpha_mat) <- list("", names(object$alpha_discount))
  
  # Print output
  cat("\n")
  cat("Call:\n")
  print(object$args1$call)
  cat("\n")
  cat("Residuals:\n")
  print(resid_sum)
  cat("\n")
  cat("Coefficients:\n")
  print(coefs)
  cat("\n")
  cat("Discount function value (alpha):\n")
  print(alpha_mat)
  cat("\n")
  if ("sigma" %in% names(object$estimates$coefficients)) {
    s <- round(bdp_regression_summary(object, "sigma", conf.int = FALSE)$estimate, 4)
    cat(paste0("Residual standard error: ", s))
    cat("\n")
  }
})
