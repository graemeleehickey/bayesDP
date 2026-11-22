# Every covariate pattern has both outcomes, within each arm and source.
logit_overlap_data <- function() {
  expand.grid(Y = 0:1, x = c(-1, 0, 1), treatment = 0:1,
              replicate = seq_len(3))[, c("Y", "x", "treatment")]
}

fit_logit_separation_case <- function(current, historical, ...) {
  bdplogit(Y ~ treatment + x, data = current, data0 = historical,
    method = "fixed", number_mcmc_alpha = 200, number_mcmc_beta = 200, ...)
}

logit_supplied_arm_priors <- list(
  prior_treatment_effect = 0, prior_control_effect = 0,
  prior_treatment_sd = 1, prior_control_sd = 1
)

test_that("zero-event and all-event groups cannot determine borrowing weights", {
  for (source in c("current", "historical")) {
    for (arm in 0:1) {
      for (outcome in 0:1) {
        d <- list(current = logit_overlap_data(), historical = logit_overlap_data())
        rows <- d[[source]]$treatment == arm
        d[[source]]$Y[rows] <- outcome
        expect_error(do.call(fit_logit_separation_case, d), paste0(
          "borrowing for the ", c("control", "treatment")[arm + 1],
          " arm: complete or quasi-complete separation detected;.*improper"))
      }
    }
  }
})

test_that("covariate separation is detected even when every group has both outcomes", {
  for (quasi in c(FALSE, TRUE)) {
    for (arm in 0:1) {
      d <- logit_overlap_data()
      rows <- d$treatment == arm
      d$Y[rows & d$x < 0] <- 0
      d$Y[rows & d$x > 0] <- 1
      if (!quasi) d$Y[rows & d$x == 0] <- 1
      expect_equal(sort(unique(d$Y[rows])), 0:1)
      expect_error(fit_logit_separation_case(d, d), paste0(
        "borrowing for the ", c("control", "treatment")[arm + 1],
        " arm: complete or quasi-complete separation detected"))
    }
  }
})

test_that("comparison models reject rank-deficient covariates", {
  d <- logit_overlap_data()
  d$duplicate <- 2 * d$x
  expect_error(bdplogit(Y ~ treatment + x + duplicate, data = d, data0 = d,
    method = "fixed"), "borrowing for the treatment arm: model matrix is rank deficient")
})

test_that("prescribed weights skip comparison data, checks, and sampling", {
  for (alpha in c(0, .35, 1)) {
    result <- discount_logit(df = stop("Comparison data should not be used"),
      discount_function = "identity", alpha_max = alpha, fix_alpha = TRUE,
      number_mcmc_alpha = stop("Comparison draws should not be requested"),
      weibull_shape = 3, weibull_scale = .135, method = "fixed")
    expect_identical(result, list(p_hat = NA_real_, alpha_discount = alpha))
  }

  current <- logit_overlap_data()
  historical <- current
  historical$Y[historical$treatment == 0] <- 0
  fit <- do.call(fit_logit_separation_case, c(
    list(current = current, historical = historical,
         fix_alpha = TRUE, alpha_max = c(.6, .25)), logit_supplied_arm_priors))
  expect_equal(fit$alpha_discount, data.frame(treatment = .6, control = .25))
  expect_equal(nrow(fit$posterior), 200)
  expect_true(all(is.finite(as.matrix(fit$posterior))))
  expect_output(summary(fit), "Coefficients:")
})

test_that("prescribed weights do not bypass historical-prior estimation checks", {
  current <- logit_overlap_data()
  historical <- current
  historical$Y[historical$treatment == 0] <- 0
  for (alpha in c(0, .5, 1)) {
    expect_error(fit_logit_separation_case(current, historical,
      fix_alpha = TRUE, alpha_max = alpha),
      "historical arm priors: complete or quasi-complete separation detected;.*infinite maximum likelihood")
  }
  partial_priors <- logit_supplied_arm_priors
  partial_priors$prior_control_sd <- NULL
  expect_error(do.call(fit_logit_separation_case, c(
    list(current = current, historical = historical, fix_alpha = TRUE), partial_priors)),
    "historical arm priors: complete or quasi-complete separation detected")

  historical <- current
  historical$x <- 0
  expect_error(fit_logit_separation_case(current, historical, fix_alpha = TRUE),
    "historical arm priors: model matrix is rank deficient")
  fit <- do.call(fit_logit_separation_case, c(
    list(current = current, historical = historical, fix_alpha = TRUE),
    logit_supplied_arm_priors))
  expect_true(all(is.finite(as.matrix(fit$posterior))))
})

test_that("current-data checks apply with estimated and prescribed weights and proper priors", {
  historical <- logit_overlap_data()
  current <- historical
  current$Y <- as.integer(current$x > 0)
  for (fixed in c(FALSE, TRUE)) {
    expect_error(do.call(fit_logit_separation_case, c(
      list(current = current, historical = historical, fix_alpha = fixed),
      logit_supplied_arm_priors)),
      "current-data logistic model: complete or quasi-complete separation detected;.*sampler requires finite")
  }

  current <- historical
  current$x <- 0
  expect_error(do.call(fit_logit_separation_case, c(
    list(current = current, historical = historical, fix_alpha = TRUE),
    logit_supplied_arm_priors)), "current-data logistic model: model matrix is rank deficient")
})

test_that("logistic model checks accept overlap across covariate units without drawing samples", {
  withr::local_seed(1701)
  rng <- .Random.seed
  for (scale in c(1e-9, 1, 1e9)) {
    d <- logit_overlap_data()
    d$x <- d$x * scale
    expect_null(bdp_check_logit_model(Y ~ treatment + x, d,
      context = "Current model", separation_reason = "Unsupported separation."))
    expect_identical(.Random.seed, rng)
  }
})

test_that("logistic prescribed-weight settings are validated", {
  d <- logit_overlap_data()
  for (alpha in list(-.1, 1.1, NA_real_, Inf, numeric(0), c(.1, .2, .3))) {
    expect_error(fit_logit_separation_case(d, d, fix_alpha = TRUE, alpha_max = alpha),
      "alpha_max must contain one or two finite values in")
  }
  for (fixed in list(NA, 1, c(TRUE, FALSE))) {
    expect_error(fit_logit_separation_case(d, d, fix_alpha = fixed), "fix_alpha must be TRUE or FALSE")
  }
})
