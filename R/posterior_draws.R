#' Convert saved posterior draws
#'
#' Use \code{posterior::as_draws_df()} or \code{posterior::as_draws_matrix()} to extract saved draws in a standard format for posterior calculations and plotting. The optional \pkg{posterior} package must be installed.
#'
#' @param x A fitted \code{bdpbinomial}, \code{bdpnormal}, \code{bdplm}, or \code{bdpsurvival} object. Fits from \code{bdplogit()} use class \code{bdplm}.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @details Conversion preserves the order and number of saved draws and does not refit, resample, or alter R's random-number state. Variable definitions are:
#' \describe{
#'   \item{Binomial and normal fits}{\code{treatment}, then \code{control} and \code{difference} when a control arm is present. These are event probabilities or population means and their paired treatment-minus-control differences, including historical-only controls and fits with \code{compare = FALSE}. The conditional mean variances stored by normal fits are not included: they are not residual-variance draws.}
#'   \item{Regression fits}{The intercept \code{(Intercept)}, treatment effect, and covariate coefficients on the original covariate scale. Logistic coefficients are on the log-odds scale. Linear fits additionally include the residual SD \code{sigma}.}
#'   \item{Survival fits}{Interval log hazards named \code{log_hazard_treatment[1]}, \code{log_hazard_treatment[2]}, and so on, with corresponding \code{log_hazard_control[j]} columns when a control arm is present. Interval order follows \code{x$args1$breaks}. One-arm fits also include \code{survival}, the survival probability at \code{x$args1$surv_time}; two-arm fits include \code{treatment}, the approximate common log hazard ratio reported by \code{tidy()}. This summary does not impose proportional hazards. Older fits without saved log hazards use the logarithms of their stored hazards.}
#' }
#' Discount weights, comparison probabilities, and separate prior or current-only draws are not included.
#'
#' The draws are represented as one sequence (\code{.chain = 1}), with consecutive \code{.iteration} and \code{.draw} indices. These indices describe retained row order, not the original iteration numbers before burn-in or thinning. Logistic fits contain a single retained MCMC sequence; the other fits contain Monte Carlo draws from the package's sampling calculations. Conversion does not create independent chains or establish convergence. In particular, splitting these rows into artificial chains is not a substitute for independent fits when assessing convergence.
#'
#' @return A \code{posterior::draws_df} or \code{posterior::draws_matrix} object, with one row per retained draw. The data-frame format includes \code{.chain}, \code{.iteration}, and \code{.draw} metadata columns.
#' @author Graeme L. Hickey
#' @seealso \code{\link{tidy}}, \code{\link{regression_extractors}}
#' @examples
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   set.seed(42)
#'   fit <- bdpbinomial(y_t = 8, N_t = 40, y_c = 5, N_c = 40,
#'     y0_t = 10, N0_t = 50, y0_c = 8, N0_c = 50,
#'     method = "fixed", number_mcmc = 100)
#'   draws <- posterior::as_draws_df(fit)
#'   mean(draws$difference > 0)
#'   posterior::as_draws_matrix(fit)
#' }
#' @name posterior_draws
NULL

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.bdpbinomial <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_matrix(as.data.frame(bdp_arm_draws(x, "posterior")))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.bdpnormal <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_matrix(as.data.frame(bdp_arm_draws(x, "posterior_mu")))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.bdplm <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_matrix(bdp_regression_draws(x, "all"))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.bdpsurvival <- function(x, ...) {
  bdp_check_dots(...)
  interval_draws <- function(arm, label) {
    draws <- arm$posterior_log_hazard
    if (is.null(draws)) draws <- log(arm$posterior_hazard)
    colnames(draws) <- paste0("log_hazard_", label, "[", seq_len(ncol(draws)), "]")
    draws
  }
  draws <- interval_draws(x$posterior_treatment, "treatment")
  if (x$args1$arm2) {
    draws <- cbind(draws, interval_draws(x$posterior_control, "control"),
      treatment = as.numeric(bdp_survival_effect(x$posterior_treatment, x$posterior_control)))
  } else {
    draws <- cbind(draws, survival = bdp_survival_probability_draws(x, x$args1$surv_time))
  }
  posterior::as_draws_matrix(draws)
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bdpbinomial <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_df(as_draws_matrix.bdpbinomial(x))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bdpnormal <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_df(as_draws_matrix.bdpnormal(x))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bdplm <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_df(as_draws_matrix.bdplm(x))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.bdpsurvival <- function(x, ...) {
  bdp_check_dots(...)
  posterior::as_draws_df(as_draws_matrix.bdpsurvival(x))
}
