#' Extract regression estimates and credible intervals
#'
#' Extract posterior summaries from linear and logistic regression fits without refitting or drawing new samples.
#'
#' @param object A fit from \code{\link{bdplm}} or \code{\link{bdplogit}}. Both use class \code{bdplm}.
#' @param parm Coefficient names or positive integer positions. Defaults to all coefficients.
#' @param level Credible interval probability, strictly between zero and one.
#' @param complete Logical, for compatibility with \code{stats::vcov()}. Both values return the same matrix because these fits do not contain aliased coefficients.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @details Coefficients are posterior means on the original covariate scale, with the intercept named \code{(Intercept)}. Logistic coefficients remain on the log-odds scale. The residual SD \code{sigma} is excluded; use \code{tidy(object, effects = "sigma")} for that parameter in a linear model.
#'
#' Despite the conventional function names, \code{confint()} returns equal-tailed posterior credible intervals, using \code{stats::quantile()} with its default interpolation. \code{vcov()} returns the sample covariance of the saved coefficient draws, including their posterior correlations. It is neither a frequentist sampling covariance nor the Monte Carlo covariance of the estimated posterior means.
#'
#' @return \code{coef()} returns a named numeric vector; \code{vcov()} returns a square matrix with matching coefficient names; \code{confint()} returns a matrix with one row per selected coefficient and two columns giving the lower and upper posterior quantiles.
#' @author Graeme L. Hickey
#' @seealso \code{\link{tidy}}, \code{\link{posterior_draws}}, \code{\link{nobs.bdpbinomial}}
#' @examples
#' set.seed(42)
#' current <- data.frame(treatment = rep(0:1, each = 20))
#' current$y <- 2 + current$treatment + rnorm(40)
#' historical <- current
#' historical$y <- 2 + historical$treatment + rnorm(40)
#' fit <- bdplm(y ~ treatment, data = current, data0 = historical,
#'   method = "fixed", number_mcmc_alpha = 100,
#'   number_mcmc_sigmagrid = 100, number_mcmc_sigma = 50,
#'   number_mcmc_beta = 200)
#' coef(fit)
#' vcov(fit)
#' confint(fit, parm = "treatment", level = 0.9)
#' @name regression_extractors
#' @importFrom stats coef confint nobs
NULL

#' @rdname regression_extractors
#' @export
coef.bdplm <- function(object, ...) {
  bdp_check_dots(...)
  colMeans(bdp_regression_draws(object))
}

#' @rdname regression_extractors
#' @export
vcov.bdplm <- function(object, complete = TRUE, ...) {
  bdp_check_dots(...)
  bdp_check_logical(complete, "complete")
  stats::cov(bdp_regression_draws(object))
}

#' @rdname regression_extractors
#' @export
confint.bdplm <- function(object, parm, level = 0.95, ...) {
  bdp_check_dots(...)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be a single finite number strictly between zero and one.", call. = FALSE)
  }
  draws <- bdp_regression_draws(object)
  if (missing(parm)) parm <- names(draws)
  valid <- length(parm) > 0L && !anyNA(parm) &&
    ((is.character(parm) && all(parm %in% names(draws))) ||
      (is.numeric(parm) && all(is.finite(parm) & parm == floor(parm) & parm >= 1 & parm <= ncol(draws))))
  if (!valid) stop("parm must contain coefficient names or positive integer positions.", call. = FALSE)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  t(vapply(draws[, parm, drop = FALSE], quantile, numeric(2), probs = probs))
}

#' Count current participants
#'
#' Return the number of current participants included in a fitted analysis.
#'
#' @param object A fitted \code{bdpbinomial}, \code{bdpnormal}, \code{bdplm}, or \code{bdpsurvival} object. Fits from \code{bdplogit()} use class \code{bdplm}.
#' @param ... Additional arguments are not supported and produce an error.
#' @details Historical participants, including historical-only controls, are excluded. These are unweighted participant counts, not effective sample sizes. Survival counts use original participant records rather than split follow-up intervals. For older survival or regression fits without the necessary records, the result is \code{NA}. Use \code{glance()} to obtain historical counts separately.
#' @return A single numeric count, matching the \code{nobs} column of \code{glance(object)}.
#' @author Graeme L. Hickey
#' @seealso \code{\link{glance}}
#' @examples
#' set.seed(42)
#' fit <- bdpbinomial(y_t = 8, N_t = 40, y0_t = 12, N0_t = 60,
#'   method = "fixed", number_mcmc = 100)
#' nobs(fit)
#' glance(fit)
#' @export
nobs.bdpbinomial <- function(object, ...) {
  bdp_check_dots(...)
  glance(object)$nobs
}

#' @rdname nobs.bdpbinomial
#' @export
nobs.bdpnormal <- function(object, ...) {
  bdp_check_dots(...)
  glance(object)$nobs
}

#' @rdname nobs.bdpbinomial
#' @export
nobs.bdplm <- function(object, ...) {
  bdp_check_dots(...)
  glance(object)$nobs
}

#' @rdname nobs.bdpbinomial
#' @export
nobs.bdpsurvival <- function(object, ...) {
  bdp_check_dots(...)
  glance(object)$nobs
}

bdp_regression_draws <- function(x, effects = "coefficients") {
  terms <- switch(effects,
    coefficients = setdiff(names(x$posterior), "sigma"),
    sigma = "sigma",
    all = names(x$posterior))
  draws <- x$posterior[, terms, drop = FALSE]
  names(draws)[names(draws) == "intercept"] <- "(Intercept)"
  draws
}
