#' Tidy posterior summaries
#'
#' Summarise saved posterior draws in a tibble, without refitting or drawing new samples. These methods also work with \code{broom::tidy()}.
#'
#' @param x A fitted \code{bdpbinomial}, \code{bdpnormal}, \code{bdplm}, or \code{bdpsurvival} object. Fits from \code{bdplogit()} use class \code{bdplm}.
#' @param conf.int Logical; include equal-tailed credible limits? Defaults to \code{FALSE}.
#' @param conf.level Credible interval probability, strictly between zero and one. Defaults to \code{0.95}.
#' @param exponentiate Logical; exponentiate logistic coefficients or the two-arm survival effect? Defaults to \code{FALSE}. Not supported for other quantities.
#' @param effects For regression fits, \code{"coefficients"} returns the intercept, treatment effect, and covariate coefficients; \code{"sigma"} returns the residual SD of a linear model. Logistic fits have no residual SD parameter.
#' @param time For a one-arm survival fit, a single finite, nonnegative time at which to summarise survival. Defaults to the fit's \code{surv_time}. Cannot be supplied for a two-arm fit.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @details
#' Binomial and normal fits return \code{treatment}, then \code{control} and \code{difference} when a control arm is present. These describe event probabilities or population means, and the treatment-minus-control difference. Estimates are posterior medians. Differences are calculated from paired draws, including when controls are historical only or the fit used \code{compare = FALSE}.
#'
#' Regression fits return posterior means on the original covariate scale, with the intercept named \code{(Intercept)}. Logistic coefficients are on the log-odds scale unless \code{exponentiate = TRUE}.
#'
#' One-arm survival fits return a \code{survival} row containing the posterior median survival probability and a \code{time} column. Two-arm fits return a \code{treatment} row for the approximate common log hazard ratio, treatment relative to control, with \code{time = NA}. Its estimate is the posterior mean of the inverse-variance-weighted interval log hazard ratio draws used by \code{summary()}; it is not a time-specific hazard ratio or an assumption of proportional hazards.
#'
#' With \code{exponentiate = TRUE}, the estimate is the exponential of the posterior mean log effect, not the posterior mean odds ratio or hazard ratio. Credible limits are exponentiated log-scale quantiles. The \code{std.error} column is the SD of the exponentiated draws, so it is on the reported ratio scale, not the log scale or a delta-method approximation. For an exponentiated logistic intercept, the estimate describes baseline odds.
#'
#' @return A tibble with one row per quantity and columns \code{term}, \code{estimate}, and \code{std.error}. Despite its conventional name, \code{std.error} is a posterior SD, not a frequentist standard error or Monte Carlo standard error. If \code{conf.int = TRUE}, \code{conf.low} and \code{conf.high} contain equal-tailed credible limits, computed using \code{stats::quantile()} with its default interpolation. Survival results also contain \code{time}. Results retain full numerical precision. No test statistics or p-values are reported.
#' @author Graeme L. Hickey
#' @seealso \code{\link{glance}}, \code{\link[=summary,bdplm-method]{summary}}
#' @examples
#' set.seed(42)
#' fit <- bdpbinomial(y_t = 8, N_t = 40, y_c = 5, N_c = 40,
#'   y0_t = 10, N0_t = 50, y0_c = 8, N0_c = 50,
#'   method = "fixed", number_mcmc = 1000)
#' tidy(fit, conf.int = TRUE)
#' glance(fit)
#' @name tidy
#' @importFrom generics tidy
#' @export
generics::tidy

#' @rdname tidy
#' @export
tidy.bdpbinomial <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  bdp_check_tidy_args(conf.int, conf.level, ...)
  bdp_arm_summary(x, "posterior", conf.int, conf.level)
}

#' @rdname tidy
#' @export
tidy.bdpnormal <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  bdp_check_tidy_args(conf.int, conf.level, ...)
  bdp_arm_summary(x, "posterior_mu", conf.int, conf.level)
}

#' @rdname tidy
#' @export
tidy.bdplm <- function(x, conf.int = FALSE, conf.level = 0.95,
                      exponentiate = FALSE, effects = c("coefficients", "sigma"), ...) {
  bdp_check_tidy_args(conf.int, conf.level, ...)
  bdp_check_logical(exponentiate, "exponentiate")
  effects <- match.arg(effects)
  if (exponentiate && (!identical(x$args1$family, "binomial") || effects != "coefficients")) {
    stop("exponentiate = TRUE is only supported for logistic coefficients.", call. = FALSE)
  }
  if (effects == "sigma" && !"sigma" %in% names(x$posterior)) {
    stop("This fit has no residual SD parameter (sigma).", call. = FALSE)
  }
  bdp_regression_summary(x, effects, conf.int, conf.level, exponentiate)
}

#' @rdname tidy
#' @export
tidy.bdpsurvival <- function(x, conf.int = FALSE, conf.level = 0.95,
                            exponentiate = FALSE, time = x$args1$surv_time, ...) {
  bdp_check_tidy_args(conf.int, conf.level, ...)
  bdp_check_logical(exponentiate, "exponentiate")
  if (x$args1$arm2) {
    if (!missing(time)) stop("time is only supported for one-arm survival fits.", call. = FALSE)
    result <- bdp_survival_effect_summary(x, conf.int, conf.level, exponentiate)
    result$time <- NA_real_
  } else {
    if (exponentiate) stop("exponentiate = TRUE is only supported for two-arm survival effects.", call. = FALSE)
    if (!is.numeric(time) || length(time) != 1L || !is.finite(time) || time < 0) {
      stop("time must be a single finite, nonnegative number.", call. = FALSE)
    }
    result <- bdp_survival_summary(x, time, conf.int, conf.level)
    result$time <- as.numeric(time)
  }
  result[, c("term", "time", setdiff(names(result), c("term", "time")))]
}

#' Glance at a Bayesian discount prior fit
#'
#' Return sample sizes, retained draw counts, and historical discount weights in a one-row tibble. These methods also work with \code{broom::glance()}.
#'
#' @param x A fitted \code{bdpbinomial}, \code{bdpnormal}, \code{bdplm}, or \code{bdpsurvival} object. Fits from \code{bdplogit()} use class \code{bdplm}.
#' @param ... Additional arguments are not supported and produce an error.
#' @return A one-row tibble containing:
#' \describe{
#'   \item{nobs}{Number of current participants included in the analysis, summed across arms.}
#'   \item{n.historical}{Number of historical participants included in the analysis, summed across arms. Historical-only controls are counted here, not in \code{nobs}. An absent source contributes zero. These are unweighted counts, not effective sample sizes.}
#'   \item{n.draws}{Number of retained outcome posterior draws, counted from the saved draws rather than the requested simulation size. This is not a simulation effective sample size.}
#'   \item{alpha.treatment, alpha.control}{Median stored historical discount weight for each arm. A prescribed or estimated fixed weight is returned unchanged. \code{NA} indicates that a weight is not defined, including historical-only controls and absent arms or historical sources.}
#'   \item{n.event, n.event.historical}{For survival fits only, numbers of observed events among current and historical participants, respectively.}
#' }
#' @details Survival sample sizes and event counts use the saved participant records, not the rows created by splitting follow-up into hazard intervals. For older survival fits without participant records, these counts are \code{NA}; refit to obtain them. The columns and their order are fixed within each fitted class. No likelihood criteria, convergence diagnostics, or posterior comparison probabilities are supplied.
#' @author Graeme L. Hickey
#' @seealso \code{\link{tidy}}
#' @examples
#' set.seed(42)
#' fit <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
#'   mu0_t = 2.8, sigma0_t = 2, N0_t = 60,
#'   method = "fixed", number_mcmc = 1000)
#' glance(fit)
#' @name glance
#' @importFrom generics glance
#' @export
generics::glance

#' @rdname glance
#' @export
glance.bdpbinomial <- function(x, ...) {
  bdp_check_dots(...)
  bdp_glance_aggregate(x, "posterior")
}

#' @rdname glance
#' @export
glance.bdpnormal <- function(x, ...) {
  bdp_check_dots(...)
  bdp_glance_aggregate(x, "posterior_mu")
}

#' @rdname glance
#' @export
glance.bdplm <- function(x, ...) {
  bdp_check_dots(...)
  data <- x$args1$data
  bdp_glance_counts(
    if (is.null(data)) NA_real_ else sum(data$historical == 0),
    if (is.null(data)) NA_real_ else sum(data$historical == 1),
    nrow(x$posterior), x$alpha_discount$treatment, x$alpha_discount$control
  )
}

#' @rdname glance
#' @export
glance.bdpsurvival <- function(x, ...) {
  bdp_check_dots(...)
  data <- x$args1$participants
  result <- bdp_glance_counts(
    if (is.null(data)) NA_real_ else sum(data$historical == 0),
    if (is.null(data)) NA_real_ else sum(data$historical == 1),
    nrow(x$posterior_treatment$posterior_hazard),
    x$posterior_treatment$alpha_discount, x$posterior_control$alpha_discount
  )
  result$n.event <- if (is.null(data)) NA_real_ else sum(data$status[data$historical == 0])
  result$n.event.historical <- if (is.null(data)) NA_real_ else sum(data$status[data$historical == 1])
  result
}

# Numerical summaries shared with the printed summary methods.
bdp_summarise_draws <- function(draws, point = median, conf.int = TRUE,
                                conf.level = 0.95, exponentiate = FALSE) {
  result <- tibble::tibble(
    term = names(draws),
    estimate = unname(vapply(draws, point, numeric(1))),
    std.error = unname(vapply(draws, function(z) sd(if (exponentiate) exp(z) else z), numeric(1)))
  )
  if (conf.int) {
    limits <- vapply(draws, quantile, numeric(2),
      probs = c((1 - conf.level) / 2, (1 + conf.level) / 2), names = FALSE)
    result$conf.low <- unname(limits[1, ])
    result$conf.high <- unname(limits[2, ])
  }
  if (exponentiate) {
    result$estimate <- exp(result$estimate)
    if (conf.int) {
      result$conf.low <- exp(result$conf.low)
      result$conf.high <- exp(result$conf.high)
    }
  }
  result
}

bdp_arm_summary <- function(x, component, conf.int = TRUE, conf.level = 0.95) {
  draws <- list(treatment = x$posterior_treatment[[component]])
  if (x$args1$arm2) {
    draws$control <- x$posterior_control[[component]]
    draws$difference <- draws$treatment - draws$control
  }
  bdp_summarise_draws(draws, conf.int = conf.int, conf.level = conf.level)
}

bdp_regression_summary <- function(x, effects = "coefficients", conf.int = TRUE,
                                   conf.level = 0.95, exponentiate = FALSE) {
  terms <- if (effects == "sigma") "sigma" else setdiff(names(x$posterior), "sigma")
  draws <- as.list(x$posterior[, terms, drop = FALSE])
  names(draws)[names(draws) == "intercept"] <- "(Intercept)"
  bdp_summarise_draws(draws, mean, conf.int, conf.level, exponentiate)
}

bdp_survival_summary <- function(x, times, conf.int = TRUE, conf.level = 0.95) {
  arm <- x$posterior_treatment
  log_hazard <- arm$posterior_log_hazard
  draws <- lapply(times, function(time) {
    if (time == 0) return(rep(1, nrow(arm$posterior_hazard)))
    if (is.null(log_hazard)) {
      lengths <- bdp_interval_lengths(time, x$args1$breaks)
      active <- lengths > 0
      return(exp(-drop(arm$posterior_hazard[, active, drop = FALSE] %*% lengths[active])))
    }
    exp(-exp(bdp_log_cumulative_hazard(log_hazard, time, x$args1$breaks)))
  })
  names(draws) <- rep("survival", length(times))
  bdp_summarise_draws(draws, conf.int = conf.int, conf.level = conf.level)
}

bdp_survival_effect_summary <- function(x, conf.int = TRUE, conf.level = 0.95,
                                       exponentiate = FALSE) {
  draws <- list(treatment = as.numeric(bdp_survival_effect(x$posterior_treatment, x$posterior_control)))
  bdp_summarise_draws(draws, mean, conf.int, conf.level, exponentiate)
}

bdp_glance_counts <- function(n, n0, draws, alpha_t, alpha_c) {
  weight <- function(z) if (length(z)) median(z) else NA_real_
  tibble::tibble(nobs = as.numeric(n), n.historical = as.numeric(n0),
    n.draws = as.numeric(draws), alpha.treatment = weight(alpha_t), alpha.control = weight(alpha_c))
}

bdp_glance_aggregate <- function(x, component) {
  bdp_glance_counts(sum(x$args1$N_t, x$args1$N_c), sum(x$args1$N0_t, x$args1$N0_c),
    length(x$posterior_treatment[[component]]),
    x$posterior_treatment$alpha_discount, x$posterior_control$alpha_discount)
}

bdp_check_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

bdp_check_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    labels <- names(dots)
    if (is.null(labels)) labels <- rep("", length(dots))
    labels[labels == ""] <- "<unnamed>"
    stop("Unsupported argument(s): ", paste(labels, collapse = ", "), ".", call. = FALSE)
  }
}

bdp_check_tidy_args <- function(conf.int, conf.level, ...) {
  bdp_check_dots(...)
  bdp_check_logical(conf.int, "conf.int")
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a single finite number strictly between zero and one.", call. = FALSE)
  }
}
