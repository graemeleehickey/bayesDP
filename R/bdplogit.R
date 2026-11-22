#' @title Bayesian Discount Prior: Two-Arm Logistic Regression
#' @description \code{bdplogit} is used to estimate the treatment effect
#'   in the presence of covariates using the logistic regression analysis
#'   implementation of the Bayesian discount prior. \code{summary} and
#'   \code{print} methods are supported. Currently, the function only supports
#'   a two-arm clinical trial where all of current treatment, current control,
#'   historical treatment, and historical control data are present
#' @param formula an object of class "formula." See "Details" for
#'   more information, including specification of treatment
#'   data indicators.
#' @param data a data frame containing the current data variables in the model.
#'   A column named \code{treatment} must be present; \code{treatment} must
#'   be binary and indicate treatment group vs. control group.
#' @param data0 a data frame containing the historical data variables in the model.
#'   The column labels of data and data0 must match.
#' @param prior_treatment_effect scalar. The historical adjusted treatment effect.
#'   If left \code{NULL}, value is estimated from the historical data.
#' @param prior_control_effect scalar. The historical adjusted control effect.
#'   If left \code{NULL}, value is estimated from the historical data.
#' @param prior_treatment_sd scalar. The standard deviation of the historical
#'   adjusted treatment effect. If left \code{NULL}, value is estimated from
#'   the historical data.
#' @param prior_control_sd scalar. The standard deviation of the historical
#'   adjusted control effect. If left \code{NULL}, value is estimated from
#'   the historical data.
#' @param prior_covariate_effect vector. The prior mean(s) of the covariate
#'   effect(s). Default value is zero. If a single value is input, the
#'   the scalar is repeated to the length of the input covariates. Otherwise,
#'   care must be taken to ensure the length of the input matches the number of
#'   covariates.
#' @param prior_covariate_sd vector. The prior standard deviation(s) of the
#'   covariate effect(s). Default value is 1e4. If a single value is input, the
#'   the scalar is repeated to the length of the input covariates. Otherwise,
#'   care must be taken to ensure the length of the input matches the number of
#'   covariates. Note that covariate effects are intentionally assigned a
#'   near-zero historical discount weight, so their priors are effectively flat:
#'   the value supplied here has negligible influence on the posterior, and the
#'   effective prior standard deviation is much larger than the nominal value
#'   (by a factor of roughly 1e6 at the default). In practice the covariate
#'   coefficients are estimated almost entirely from the current data, and the
#'   posterior is insensitive to this argument across many orders of magnitude.
#' @param discount_function character. Specify the discount function to use.
#'   Currently supports \code{weibull}, \code{scaledweibull}, and
#'   \code{identity}. The discount function \code{scaledweibull} scales
#'   the output of the Weibull CDF to have a max value of 1. The \code{identity}
#'   discount function uses the posterior probability directly as the discount
#'   weight. Default value is "\code{identity}".
#' @param alpha_max scalar in [0, 1]. Maximum weight the discount function can apply.
#'   Default is 1. Users may specify a vector of two values where the first
#'   value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#'   If TRUE, comparison models are not fitted and the prescribed weights are
#'   used directly. Historical-prior estimation and current-data model checks
#'   still apply.
#' @param number_mcmc_alpha scalar. Number of Monte Carlo
#'   simulations for estimating the historical data weight. Default is 5000.
#' @param number_mcmc_beta scalar. Number of Monte Carlo simulations for
#'   estimating beta, the vector of regression coefficients. Default is 10000.
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. Users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group. Not used when \code{discount_function} = "identity".
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. Users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group.
#'   Not used when \code{discount_function} = "identity".
#' @param method character. Set to \code{"fixed"}, the only supported method.
#'   It estimates each historical weight once and holds it fixed throughout
#'   the analysis, unless \code{fix_alpha = TRUE} prescribes the weights.
#'   The legacy default \code{"mc"} is not supported and raises an error.
#' @details \code{bdplogit} uses a two-stage approach for determining the
#'   strength of historical data in estimation of an adjusted mean or covariate
#'   effect. In the first stage, a \emph{discount function} is used that
#'   defines the maximum strength of the
#'   historical data and discounts based on disagreement with the current data.
#'   Disagreement between current and historical data is assessed by fitting
#'   a logistic regression separately within each arm, with a historical
#'   indicator and the covariates, using flat coefficient priors. If \eqn{q}
#'   is the posterior probability that the historical-indicator coefficient is
#'   positive, the comparison is \eqn{p = 2\min(q, 1-q)}. The discount function
#'   transforms \eqn{p} into the historical weight \code{alpha}. This is not
#'   a posterior probability of equality or equivalence.
#'
#'   In the second stage, posterior estimation is performed where the discount
#'   function parameter, \code{alpha}, is used to weight the historical data
#'   effects.
#'
#'   The formula must include an intercept (i.e., do not use \code{-1} in
#'   the formula) and both data and data0 must be present.
#'   The column names of data and data0 must match. See \code{examples} below for
#'   example usage.
#'
#'   Formula offsets are not supported. Calls to \code{offset()} or
#'   \code{stats::offset()} cause an error before fitting, including offsets
#'   that are zero or cancel each other. Prescribing weights, disabling
#'   borrowing, or supplying historical priors does not bypass this restriction.
#'
#'   Internally, the model is parameterized without a fixed intercept, using
#'   separate treatment and control (log-odds) means rather than an intercept
#'   plus a treatment effect. In this parameterization the arm means are the
#'   fitted values at covariate = 0, so when covariates are not centered the
#'   treatment- and control-mean estimators become strongly correlated and
#'   their standard errors are extrapolation errors at covariate = 0. This
#'   distorts the discount-prior construction and the historical borrowing. To
#'   guard against this, \code{bdplogit} automatically mean-centers each
#'   covariate on its pooled (current plus historical) mean before fitting, and
#'   back-transforms the reported intercept onto the original covariate scale.
#'   Specifically, if \eqn{\bar{x}_j} is the pooled mean of covariate \eqn{j}
#'   and \eqn{\beta_j} its estimated effect, the reported intercept is
#'   \eqn{\beta_0 = \beta_0^{c} - \sum_j \beta_j \bar{x}_j}, where
#'   \eqn{\beta_0^{c}} is the control-arm log-odds estimated on the centered
#'   scale. The treatment effect and covariate slopes are unchanged by
#'   centering; the reported \code{intercept} remains the control-arm log-odds
#'   at covariate = 0 on the original scale. As a result, estimates are
#'   invariant to a location shift of any covariate, and users do not need to
#'   center covariates themselves.
#'
#'   Posterior sampling uses \code{MCMCpack::MCMClogit()}. Before fitting,
#'   model matrices are checked for rank deficiency, and complete or
#'   quasi-complete separation is detected using \code{detectseparation}.
#'   These checks include separation caused by covariates, as well as groups
#'   with no events or with all events. An error identifies the affected stage:
#'   \itemize{
#'     \item Each arm's comparison model must have a proper posterior under
#'       its flat coefficient priors. Separation or rank deficiency prevents
#'       estimation of that arm's borrowing weight.
#'     \item When any historical arm-prior mean or SD is estimated from the
#'       historical regression, that regression must have a full-rank design
#'       and no separation. Numerically large but finite estimates from a
#'       separated regression are not used as historical priors.
#'     \item The current-data regression must also have a full-rank design
#'       and no separation because the sampler uses an unpenalized logistic
#'       regression to construct starting values and proposals. This is a
#'       limitation of the supported sampling procedure, not a claim that
#'       separation makes a posterior with proper coefficient priors improper.
#'   }
#'   With \code{fix_alpha = TRUE}, the comparison models and their checks are
#'   skipped; their internal \code{p_hat} values are \code{NA_real_}.
#'   Historical regression checks are skipped only when all four arm-prior
#'   means and SDs are supplied. Current-data checks always apply. Prescribing
#'   weights alone therefore does not make separated datasets supported.
#'
#' @return \code{bdplogit} returns an object of class "bdplm", reusing its
#'   print and summary methods. Residual summaries use observed binary outcomes
#'   minus fitted probabilities from the posterior mean coefficients. The
#'   centering means are saved in \code{args1$covariate_means}; older saved fits
#'   without these means must be refitted before summarising residuals.
#'
#' The returned object is a list containing at least
#' the following components:
#' \describe{
#'   \item{\code{posterior}}{
#'     data frame. Posterior draws of the intercept, treatment effect, and
#'     covariate effects on the log-odds scale.}
#'   \item{\code{alpha_discount}}{
#'     data frame. One row of estimated or prescribed historical weights,
#'     with columns \code{treatment} and \code{control}. These weights are
#'     not posterior probabilities of equality.}
#'   \item{\code{estimates}}{
#'     list. Posterior means and posterior SDs of the intercept, treatment
#'     effect, and covariate effects.}
#' }
#'
#' @examples
#' # Set sample sizes
#' n_t <- 30 # Current treatment sample size
#' n_c <- 30 # Current control sample size
#' n_t0 <- 80 # Historical treatment sample size
#' n_c0 <- 80 # Historical control sample size
#'
#' # Treatment group vectors for current and historical data
#' treatment <- c(rep(1, n_t), rep(0, n_c))
#' treatment0 <- c(rep(1, n_t0), rep(0, n_c0))
#'
#' # Simulate a covariate effect for current and historical data
#' x <- rnorm(n_t + n_c, 1, 1)
#' x0 <- rnorm(n_t0 + n_c0, 1, 1)
#'
#' # Simulate a binary outcome on the logit scale:
#' # - Intercept of -1 for current and historical data
#' # - Treatment (log-odds) effect of 1.0 for current data
#' # - Treatment (log-odds) effect of 0.9 for historical data
#' # - Covariate effect of 0.3 for current and historical data
#' Y <- rbinom(n_t + n_c, 1, plogis(-1 + 1.0 * treatment + 0.3 * x))
#' Y0 <- rbinom(n_t0 + n_c0, 1, plogis(-1 + 0.9 * treatment0 + 0.3 * x0))
#'
#' # Place data into separate treatment and control data frames and
#' # assign historical = 0 (current) or historical = 1 (historical)
#' df_ <- data.frame(Y = Y, treatment = treatment, x = x)
#' df0 <- data.frame(Y = Y0, treatment = treatment0, x = x0)
#'
#' # Fit model using default settings
#' fit <- bdplogit(
#'   formula = Y ~ treatment + x, data = df_, data0 = df0,
#'   method = "fixed"
#' )
#'
#' # Look at estimates and discount weight
#' summary(fit)
#' print(fit)
#' @rdname bdplogit
#' @import methods
#' @importFrom stats density is.empty.model median model.offset model.response
#'   pweibull pnorm quantile rbeta rgamma rnorm var vcov contrasts<- dt gaussian
#'   lm.fit model.frame model.matrix.default offset terms terms.formula
#'   coefficients lm qgamma runif glm binomial
#' @importFrom MCMCpack MCMClogit
#' @aliases bdplogit-method
#' @aliases bdplogit,ANY-method
#' @useDynLib bayesDP, .registration = TRUE
#' @export bdplogit
bdplogit <- setClass("bdplogit", slots = c(
  posterior_treatment = "list",
  posterior_control = "list",
  args1 = "list"
))
setGeneric(
  "bdplogit",
  function(formula = formula,
           data = data,
           data0 = NULL,
           prior_treatment_effect = NULL, # 0,
           prior_control_effect = NULL, # 0,
           prior_treatment_sd = NULL, # 1e4,
           prior_control_sd = NULL, # 1e4,
           prior_covariate_effect = 0,
           prior_covariate_sd = 1e4,
           number_mcmc_alpha = 5000,
           number_mcmc_beta = 10000,
           discount_function = "identity",
           alpha_max = 1,
           fix_alpha = FALSE,
           weibull_scale = 0.135,
           weibull_shape = 3,
           method = "mc") {
    standardGeneric("bdplogit")
  }
)

setMethod(
  "bdplogit",
  signature(),
  function(formula = formula,
           data = data,
           data0 = NULL,
           prior_treatment_effect = NULL, # 0,
           prior_control_effect = NULL, # 0,
           prior_treatment_sd = NULL, # 1e4,
           prior_control_sd = NULL, # 1e4,
           prior_covariate_effect = 0,
           prior_covariate_sd = 1e4,
           number_mcmc_alpha = 5000,
           number_mcmc_beta = 10000,
           discount_function = "identity",
           alpha_max = 1,
           fix_alpha = FALSE,
           weibull_scale = 0.135,
           weibull_shape = 3,
           method = "mc") {

    ### Check validity of data input
    call <- match.call()
    if (missing(data)) {
      stop("Current data not input correctly.")
    }

    if (is.null(data0)) {
      stop("Historical data not input correctly.")
    }

    if (method == "mc") {
      stop("Method 'mc' not currently supported.")
    }
    if (!is.logical(fix_alpha) || length(fix_alpha) != 1L || is.na(fix_alpha)) {
      stop("fix_alpha must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.numeric(alpha_max) || !length(alpha_max) %in% c(1L, 2L) ||
        any(!is.finite(alpha_max)) || any(alpha_max < 0 | alpha_max > 1)) {
      stop("alpha_max must contain one or two finite values in [0, 1].", call. = FALSE)
    }

    ### Parse current data
    bdp_check_no_offsets(formula, data, "bdplogit")
    mf <- mf0 <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action <- NULL
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) {
        names(Y) <- nm
      }
    }

    X <- if (!is.empty.model(mt)) {
      model.matrixBayes(
        object = mt, data = data, contrasts.arg = NULL,
        keep.order = TRUE, drop.baseline = TRUE
      )
    } else {
      matrix(, NROW(Y), 0L)
    }

    ### Alter intercept column name, if present
    imatch <- match("(Intercept)", colnames(X))
    if (!is.na(imatch)) colnames(X)[imatch] <- "intercept"

    # Create indicator of whether each column is present
    cmatch <- match(c("intercept", "treatment"), colnames(X))
    cmatch <- !is.na(cmatch) ### == TRUE == present

    if (!all(cmatch)) {
      stop("Current data is input incorrectly. Intercept and treatment must be present.")
    }

    ### Count levels of treatment data and sure 0 and 1 are present
    trt_levels <- levels(as.factor(X[, "treatment"]))

    if (!(all(trt_levels %in% c(0, 1)))) {
      stop("Treatment input has wrong levels. Values should be 0 and 1.")
    }

    ### Parse historical data
    if (missing(data0)) {
      stop("Historical data is required.")
    } else {
      m00 <- match("data0", names(mf0), 0L)
      m01 <- match("data", names(mf0), 0L)
      mf0[m01] <- mf0[m00]
      m0 <- match(c("formula", "data"), names(mf0), 0L)
      mf0 <- mf0[c(1L, m0)]
      mf0$na.action <- NULL
      mf0[[1L]] <- quote(stats::model.frame)
      mf0 <- eval(mf0, parent.frame())
      mt0 <- attr(mf0, "terms")
      Y0 <- model.response(mf0, "any")
      if (length(dim(Y0)) == 1L) {
        nm0 <- rownames(Y0)
        dim(Y0) <- NULL
        if (!is.null(nm0)) {
          names(Y0) <- nm0
        }
      }

      X0 <- if (!is.empty.model(mt0)) {
        model.matrixBayes(
          object = mt0, data = data0, contrasts.arg = NULL,
          keep.order = TRUE, drop.baseline = TRUE
        )
      } else {
        matrix(, NROW(Y0), 0L)
      }

      ### Alter intercept column name, if present
      imatch <- match("(Intercept)", colnames(X0))
      if (!is.na(imatch)) colnames(X0)[imatch] <- "intercept"

      # Create indicator of whether each column is present
      cmatch <- match(c("intercept", "treatment"), colnames(X0))
      cmatch <- !is.na(cmatch) ### == TRUE == present

      if (!all(cmatch)) {
        stop("Historical data is input incorrectly. Intercept and treatment must be present.")
      }

      cnames <- colnames(X)
      cnames0 <- colnames(X0)
      if (!all(cnames == cnames0)) {
        stop("Column names in the historical data must match the current data.")
      }

      ### Count levels of treatment data and ensure 0 and 1 are present
      trt_levels0 <- levels(as.factor(X0[, "treatment"]))

      # Check that the levels are 0 and/or 1
      if (!(all(trt_levels0 %in% c(0, 1)))) {
        stop("Treatment input has wrong levels. Values should be 0 and 1.")
      }
    }

    ### If method == "mc", fix_alpha must be FALSE
    if (method == "mc" & fix_alpha) {
      stop("mc method not possible with fix_alpha == TRUE. Set fix_alpha == FALSE.")
    }

    # Check that discount_function is input correctly
    all_functions <- c("weibull", "scaledweibull", "identity")
    function_match <- match(discount_function, all_functions)
    if (is.na(function_match)) {
      stop("discount_function input incorrectly.")
    }

    historical <- NULL
    treatment <- NULL
    intercept <- NULL

    ##############################################################################
    # Quick check, if alpha_max, weibull_scale, or weibull_shape have length 1,
    # repeat input twice
    ##############################################################################

    if (length(alpha_max) == 1) {
      alpha_max <- rep(alpha_max, 2)
    }

    if (length(weibull_scale) == 1) {
      weibull_scale <- rep(weibull_scale, 2)
    }

    if (length(weibull_shape) == 1) {
      weibull_shape <- rep(weibull_shape, 2)
    }

    ##############################################################################
    # Format input data
    ##############################################################################

    ### Create a master dataframe
    df <- data.frame(Y = Y, data.frame(X), historical = 0)
    df0 <- data.frame(Y = Y0, data.frame(X0), historical = 1)

    df <- rbind(df, df0)

    ### Count number of covariates
    names_df <- names(df)
    covs_df <- names_df[!(names_df %in% c("Y", "intercept", "treatment", "historical"))]
    n_covs <- length(covs_df)

    ##############################################################################
    # Mean-center covariates
    #
    # The model is parameterized without a fixed intercept, using separate
    # treatment and control (log-odds) means (i.e., Y ~ -1 + treatment +
    # control + covs). In this parameterization the arm means are the fitted
    # values at covariate = 0. When covariates are not centered these arm-mean
    # estimators become strongly correlated and their standard errors are
    # extrapolation errors at covariate = 0, which breaks the (diagonal)
    # discount-prior construction and corrupts the historical borrowing.
    # Centering each covariate on its pooled (current + historical) mean
    # restores the intended behavior. The centering is reversed on the
    # posterior intercept below so that reported quantities remain on the
    # original covariate scale.
    ##############################################################################

    covariate_means <- vapply(covs_df, function(nm) mean(df[[nm]]), numeric(1))

    for (nm in covs_df) {
      df[[nm]] <- df[[nm]] - covariate_means[[nm]]
    }

    ### Split data into separate treatment & control dataframes
    df_t <- subset(df, treatment == 1)
    df_c <- subset(df, treatment == 0)

    ### Also create current data dataframe
    df_current <- subset(df, historical == 0, select = -intercept)

    ##############################################################################
    # Estimate discount weights for each of the treatment and control arms
    ##############################################################################

    discount_treatment <- discount_logit(
      df = df_t,
      discount_function = discount_function,
      alpha_max = alpha_max[1],
      fix_alpha = fix_alpha,
      number_mcmc_alpha = number_mcmc_alpha,
      weibull_shape = weibull_shape[1],
      weibull_scale = weibull_scale[1],
      method = method,
      arm = "treatment"
    )

    discount_control <- discount_logit(
      df = df_c,
      discount_function = discount_function,
      alpha_max = alpha_max[2],
      fix_alpha = fix_alpha,
      number_mcmc_alpha = number_mcmc_alpha,
      weibull_shape = weibull_shape[2],
      weibull_scale = weibull_scale[2],
      method = method,
      arm = "control"
    )

    ##############################################################################
    # Estimate historical adjusted treatment and control effects to use as the
    # priors for the current data
    ##############################################################################

    if (is.null(prior_treatment_effect) | is.null(prior_control_effect) |
        is.null(prior_treatment_sd) | is.null(prior_control_sd)) {
      ### Use the centered historical data so the prior arm means are estimated
      ### on the same (centered) covariate scale as the current-data design.
      df0_centered <- subset(df, historical == 1)
      df0_centered$control <- 1 - df0_centered$treatment

      cnames0 <- names(df0_centered)
      cnames0 <- cnames0[!(cnames0 %in% c("Y", "intercept", "historical", "treatment", "control"))]
      cnames0 <- c("treatment", "control", cnames0)
      f0 <- paste0("Y~ -1+", paste0(cnames0, collapse = "+"))

      bdp_check_logit_model(f0, df0_centered,
        context = "Cannot estimate historical arm priors",
        separation_reason = "the historical regression has infinite maximum likelihood estimates.")
      fit_0 <- glm(f0, data = df0_centered, family = binomial)
      summ_0 <- summary(fit_0)

      if (is.null(prior_treatment_effect)) prior_treatment_effect <- summ_0$coef[1, 1]
      if (is.null(prior_control_effect)) prior_control_effect <- summ_0$coef[2, 1]

      if (is.null(prior_treatment_sd)) prior_treatment_sd <- summ_0$coef[1, 2]
      if (is.null(prior_control_sd)) prior_control_sd <- summ_0$coef[2, 2]
    }

    ### Prior covariate effects
    if (length(prior_covariate_effect) == 1) {
      prior_covariate_effect <- rep(prior_covariate_effect, n_covs)
    }

    if (length(prior_covariate_sd) == 1) {
      prior_covariate_sd <- rep(prior_covariate_sd, n_covs)
    }

    ##############################################################################
    # Estimate augmented treatment effect
    ##############################################################################

    ### Compute prior terms
    ### - Covariance/variance needs to be parameterized as precision
    tau2 <- c(1 / prior_treatment_sd, 1 / prior_control_sd, 1 / prior_covariate_sd)^2
    mu0 <- c(prior_treatment_effect, prior_control_effect, prior_covariate_effect)

    ### Calculate constants from current data
    df_current$control <- 1 - df_current$treatment

    # Create analysis formula. The data frame must retain the response Y so
    # that MCMClogit can resolve it from the formula.
    predictors <- c("treatment", "control", covs_df)
    df_analysis <- df_current[, c("Y", predictors)]
    f <- paste0("Y ~ -1 +", paste0(predictors, collapse = "+"))

    bdp_check_logit_model(f, df_analysis,
      context = "Cannot fit the current-data logistic model",
      separation_reason = "the sampler requires finite unpenalized logistic estimates for its starting values and proposals.")

    ### Estimation scheme differs conditional on discounting method
    if (method == "fixed") {

      ### Extract alpha0, append "zero" weight for the covariate effect(s)
      alpha0 <- c(
        discount_treatment$alpha_discount + 1e-12,
        discount_control$alpha_discount + 1e-12,
        rep(1e-12, n_covs)
      )

      ### Create the discounted prior precision matrix. tau2 is already a
      ### precision (1 / sd^2), so the discount weight multiplies it.
      B0 <- diag(alpha0 * tau2)

      ### Get Bayesian fit
      fit <- MCMCpack::MCMClogit(f,
                                 data = df_analysis,
                                 mcmc = number_mcmc_beta,
                                 b0 = mu0,
                                 B0 = B0
      )

      beta_samples <- data.frame(fit)
    } else if (method == "mc") {
      stop("Method 'mc' not currently supported.")
    }

    ### Estimate posterior of intercept and treatment effect
    ### The control coefficient is the control-arm (log-odds) mean at centered
    ### covariate = 0, i.e. at the original covariate means. Back-transform to
    ### the original scale so the reported intercept is the control log-odds at
    ### covariate = 0. The treatment effect (difference of arm means) and the
    ### covariate slopes are unaffected by centering.
    beta_samples$intercept <- beta_samples$control
    if (n_covs > 0) {
      covariate_shift <- as.matrix(beta_samples[, covs_df, drop = FALSE]) %*%
        covariate_means[covs_df]
      beta_samples$intercept <- beta_samples$intercept - as.vector(covariate_shift)
    }
    beta_samples$treatment <- beta_samples$treatment - beta_samples$control

    ### Format posterior table
    m <- ncol(beta_samples)
    beta_samples <- beta_samples[, c(m, 1, 3:(m - 1))]

    ### Format alpha_discount values
    alpha_discount <- data.frame(
      treatment = discount_treatment$alpha_discount,
      control = discount_control$alpha_discount
    )

    ### Format estimates
    estimates <- list()
    estimates$coefficients <- data.frame(t(colMeans(beta_samples)))
    estimates$coefficients <- estimates$coefficients[c("intercept", "treatment", covs_df)]
    estimates$se <- data.frame(t(apply(beta_samples, 2, sd)))
    estimates$se <- estimates$se[c("intercept", "treatment", covs_df)]

    ### Format input arguments
    args1 <- list(
      call = call,
      data = df,
      covariate_means = covariate_means,
      family = "binomial"
    )

    ### Format output
    me <- list(
      posterior = beta_samples,
      alpha_discount = alpha_discount,
      estimates = estimates,
      args1 = args1
    )

    ### bdplogit reuses the bdplm S4 print/summary methods, since the output
    ### structure is identical. The object is therefore classed as "bdplm".
    class(me) <- "bdplm"
    return(me)
  }
)


################################################################################
# Logistic Regression discount weight estimation
# 1) Estimate the discount function
#    - Compare the historical-indicator coefficient with zero
#    - Estimate alpha from that directional posterior comparison
################################################################################

discount_logit <- function(df, discount_function, alpha_max, fix_alpha,
                           number_mcmc_alpha,
                           weibull_shape, weibull_scale,
                           method, arm = NULL) {

  if (fix_alpha) {
    return(list(p_hat = NA_real_, alpha_discount = alpha_max))
  }

  # Create formula
  cnames <- names(df)
  cnames <- cnames[!(cnames %in% c("Y", "intercept", "historical", "treatment"))]
  cnames <- c("historical", cnames)
  f <- paste0("Y~", paste0(cnames, collapse = "+"))

  context <- if (is.null(arm)) {
    "Cannot estimate historical borrowing"
  } else {
    paste0("Cannot estimate historical borrowing for the ", arm, " arm")
  }
  bdp_check_logit_model(f, df, context,
    separation_reason = "the comparison posterior under flat coefficient priors is improper.")

  ### Get Bayesian fit
  fit <- MCMCpack::MCMClogit(f, data = df, mcmc = number_mcmc_alpha, B0 = 0)
  posterior <- data.frame(fit)

  ### Monte Carlo simulations of historical effect
  beta <- posterior$historical

  ### Compute posterior comparison
  if (method == "fixed") {
    p_hat <- mean(beta > 0)
    p_hat <- 2 * ifelse(p_hat > 0.5, 1 - p_hat, p_hat)
  } else if (method == "mc") {
    stop("Method 'mc' not currently supported.")
    # Z     <- abs(beta)/sigma2_beta
    # p_hat <- 2*(1-pnorm(Z))
  }

  ### Compute alpha, the discount weight
  if (discount_function == "weibull") {
    alpha_discount <- pweibull(p_hat,
                               shape = weibull_shape,
                               scale = weibull_scale
    ) * alpha_max
  } else if (discount_function == "scaledweibull") {
    max_p <- pweibull(1, shape = weibull_shape, scale = weibull_scale)

    alpha_discount <- pweibull(p_hat,
                               shape = weibull_shape,
                               scale = weibull_scale
    ) * alpha_max / max_p
  } else if (discount_function == "identity") {
    alpha_discount <- p_hat * alpha_max
  }

  res <- list(
    p_hat = p_hat,
    alpha_discount = alpha_discount
  )
  res
}

# Check the exact model used at each logistic fitting stage before estimation.
bdp_check_logit_model <- function(formula, data, context, separation_reason) {
  mf <- stats::model.frame(stats::as.formula(formula), data = data,
                          na.action = stats::na.fail)
  X <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)
  if (!nrow(X) || any(!is.finite(X)) ||
      !(is.numeric(y) || is.logical(y)) || any(!y %in% c(0, 1))) {
    stop(context, ": finite predictors and binary outcomes coded 0/1 are required.", call. = FALSE)
  }

  # Rescale columns for the numerical checks without changing the fitted model.
  column_scale <- apply(abs(X), 2, max)
  column_scale[column_scale == 0] <- 1
  X <- sweep(X, 2, column_scale, "/")
  if (qr(X)$rank < ncol(X)) {
    stop(context, ": model matrix is rank deficient; check treatment groups and redundant covariates.", call. = FALSE)
  }

  check <- detectseparation::detect_separation(x = X, y = as.numeric(y),
    family = stats::binomial(), intercept = attr(attr(mf, "terms"), "intercept") == 1,
    control = detectseparation::detect_separation_control(
      implementation = "lpSolveAPI", purpose = "test"))
  if (!is.logical(check$outcome) || length(check$outcome) != 1L || is.na(check$outcome)) {
    stop(context, ": separation could not be determined.", call. = FALSE)
  }
  if (check$outcome) {
    stop(context, ": complete or quasi-complete separation detected; ",
         separation_reason, call. = FALSE)
  }
  invisible(NULL)
}
