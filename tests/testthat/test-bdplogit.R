## Tests for bdplogit() input validation and the fixed-effect fit path.

make_logit_data <- function() {
  set.seed(88)
  n_t <- 30
  n_c <- 30
  n_t0 <- 80
  n_c0 <- 80
  treatment <- c(rep(1, n_t), rep(0, n_c))
  treatment0 <- c(rep(1, n_t0), rep(0, n_c0))
  x <- rnorm(n_t + n_c, 1, 1)
  x0 <- rnorm(n_t0 + n_c0, 1, 1)
  Yb <- rbinom(n_t + n_c, 1, plogis(-1 + treatment + 0.3 * x))
  Yb0 <- rbinom(n_t0 + n_c0, 1, plogis(-1 + 0.9 * treatment0 + 0.3 * x0))
  list(
    data = data.frame(Y = Yb, treatment = treatment, x = x),
    data0 = data.frame(Y = Yb0, treatment = treatment0, x = x0)
  )
}

test_that("bdplogit errors when current data are missing", {
  d <- make_logit_data()
  expect_error(
    bdplogit(Y ~ treatment + x, data0 = d$data0),
    "Current data not input correctly"
  )
})

test_that("bdplogit errors when historical data are missing", {
  d <- make_logit_data()
  expect_error(
    bdplogit(Y ~ treatment + x, data = d$data),
    "Historical data not input correctly"
  )
})

test_that("bdplogit errors for the unsupported mc method", {
  d <- make_logit_data()
  expect_error(
    bdplogit(Y ~ treatment + x, data = d$data, data0 = d$data0, method = "mc"),
    "Method 'mc' not currently supported"
  )
})

test_that("bdplogit errors for an unknown discount function", {
  d <- make_logit_data()
  expect_error(
    bdplogit(Y ~ treatment + x,
      data = d$data, data0 = d$data0,
      discount_function = "not_a_function", method = "fixed"
    ),
    "discount_function input incorrectly"
  )
})

test_that("bdplogit errors when treatment levels are not 0/1", {
  d <- make_logit_data()
  bad <- d$data
  bad$treatment <- bad$treatment + 1 # levels become 1/2
  expect_error(
    bdplogit(Y ~ treatment + x, data = bad, data0 = d$data0, method = "fixed"),
    "Treatment input has wrong levels"
  )
})

test_that("bdplogit fixed-effect fit recovers the treatment effect", {
  d <- make_logit_data()
  set.seed(204)
  fit <- bdplogit(Y ~ treatment + x,
    data = d$data, data0 = d$data0, method = "fixed"
  )

  # The fitted object reuses the bdplm class for its print/summary methods.
  expect_s3_class(fit, "bdplm")
  expect_named(
    fit$estimates$coefficients,
    c("intercept", "treatment", "x")
  )

  # Treatment effect on the log-odds scale should be near the simulated 1.0.
  expect_equal(fit$estimates$coefficients$treatment, 1.0, tolerance = 0.6)

  # Discount weights are valid probabilities.
  alpha <- unlist(fit$alpha_discount)
  expect_true(all(alpha >= 0 & alpha <= 1))

  # print/summary (inherited from bdplm) run without error.
  expect_output(summary(fit), "Coefficients:")
  expect_output(print(fit), "Coefficients:")
})

test_that("bdplogit honours user-specified priors", {
  d <- make_logit_data()
  set.seed(77)
  fit <- bdplogit(Y ~ treatment + x,
    data = d$data, data0 = d$data0, method = "fixed",
    prior_treatment_effect = 0, prior_control_effect = 0,
    prior_treatment_sd = 1, prior_control_sd = 1
  )
  expect_s3_class(fit, "bdplm")
  expect_true(is.finite(fit$estimates$coefficients$treatment))
})

test_that("bdplogit treatment effect is invariant to covariate location shifts (#1)", {
  # bdplogit uses an intercept-free parameterization (separate treatment and
  # control log-odds means), so it now mean-centers covariates internally. A
  # pure location shift of a covariate only relabels where covariate = 0 sits
  # and must not change the treatment (log-odds) effect, its posterior SD, or
  # the covariate slope. The intercept (control log-odds at covariate = 0) is
  # expected to change and is checked to follow the exact back-transformation.
  set.seed(6126)
  n <- 120
  n0 <- 300
  trt <- rep(0:1, each = n / 2)
  trt0 <- rep(0:1, each = n0 / 2)
  x_base <- rnorm(n, 0, 2)
  x0_base <- rnorm(n0, 0, 2)
  Yb <- rbinom(n, 1, plogis(-0.5 + 1.2 * trt + 0.8 * x_base))
  Yb0 <- rbinom(n0, 1, plogis(-0.5 + 1.2 * trt0 + 0.8 * x0_base))

  fit_at_shift <- function(shift) {
    cur <- data.frame(Y = Yb, treatment = trt, x = x_base + shift)
    hist <- data.frame(Y = Yb0, treatment = trt0, x = x0_base + shift)
    set.seed(4471)
    bdplogit(Y ~ treatment + x, data = cur, data0 = hist, method = "fixed")
  }

  fit0 <- fit_at_shift(0)
  fit_shift <- fit_at_shift(100)

  # Treatment effect, its posterior SD, and the slope are invariant.
  expect_equal(
    mean(fit_shift$posterior$treatment),
    mean(fit0$posterior$treatment),
    tolerance = 1e-8
  )
  expect_equal(
    sd(fit_shift$posterior$treatment),
    sd(fit0$posterior$treatment),
    tolerance = 1e-8
  )
  expect_equal(
    mean(fit_shift$posterior$x),
    mean(fit0$posterior$x),
    tolerance = 1e-8
  )

  # The intercept is the control log-odds at covariate = 0, so shifting x by
  # 100 moves it by exactly -100 * slope.
  expect_equal(
    mean(fit_shift$posterior$intercept),
    mean(fit0$posterior$intercept) - 100 * mean(fit0$posterior$x),
    tolerance = 1e-6
  )
})
