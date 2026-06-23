## Tests for factor-covariate handling in model.matrixBayes() via bdplm()
## and bdplogit(). These exercise the contrast-generation path used for
## unordered factors.

make_factor_lm_data <- function() {
  set.seed(1847)
  n_t <- 60
  n_c <- 60
  n_t0 <- 120
  n_c0 <- 120
  treatment <- c(rep(1, n_t), rep(0, n_c))
  treatment0 <- c(rep(1, n_t0), rep(0, n_c0))
  site <- factor(rep(c("A", "B", "C"), length.out = n_t + n_c),
    levels = c("A", "B", "C")
  )
  site0 <- factor(rep(c("A", "B", "C"), length.out = n_t0 + n_c0),
    levels = c("A", "B", "C")
  )
  x <- rnorm(n_t + n_c)
  x0 <- rnorm(n_t0 + n_c0)
  site_eff <- c(A = 0, B = 2, C = -1)
  Y <- 5 + 4 * treatment + 1.5 * x + site_eff[as.character(site)] +
    rnorm(n_t + n_c, 0, 2)
  Y0 <- 5 + 3.5 * treatment0 + 1.5 * x0 + site_eff[as.character(site0)] +
    rnorm(n_t0 + n_c0, 0, 2)

  list(
    data = data.frame(Y = Y, treatment = treatment, x = x, site = site),
    data0 = data.frame(Y = Y0, treatment = treatment0, x = x0, site = site0)
  )
}

make_factor_logit_data <- function() {
  set.seed(2486)
  n_t <- 80
  n_c <- 80
  n_t0 <- 160
  n_c0 <- 160
  treatment <- c(rep(1, n_t), rep(0, n_c))
  treatment0 <- c(rep(1, n_t0), rep(0, n_c0))
  site <- factor(rep(c("A", "B", "C"), length.out = n_t + n_c),
    levels = c("A", "B", "C")
  )
  site0 <- factor(rep(c("A", "B", "C"), length.out = n_t0 + n_c0),
    levels = c("A", "B", "C")
  )
  x <- rnorm(n_t + n_c)
  x0 <- rnorm(n_t0 + n_c0)
  site_eff <- c(A = 0, B = 0.4, C = -0.3)
  lp <- -1 + 0.9 * treatment + 0.3 * x + site_eff[as.character(site)]
  lp0 <- -1 + 0.8 * treatment0 + 0.3 * x0 + site_eff[as.character(site0)]
  Y <- rbinom(n_t + n_c, 1, plogis(lp))
  Y0 <- rbinom(n_t0 + n_c0, 1, plogis(lp0))

  list(
    data = data.frame(Y = Y, treatment = treatment, x = x, site = site),
    data0 = data.frame(Y = Y0, treatment = treatment0, x = x0, site = site0)
  )
}

test_that("bdplm supports unordered factor covariates", {
  d <- make_factor_lm_data()
  set.seed(1921)
  fit <- bdplm(Y ~ treatment + x + site,
    data = d$data, data0 = d$data0,
    method = "fixed", number_mcmc_beta = 2000
  )

  expect_named(
    fit$estimates$coefficients,
    c("intercept", "treatment", "x", "siteB", "siteC", "sigma")
  )
  expect_true(all(is.finite(unlist(fit$estimates$coefficients))))
  expect_true(all(c("siteB", "siteC") %in% names(fit$posterior)))

  # The simulated treatment effect is 4; this tolerance allows Monte Carlo and
  # sample variation while still checking the factor model was fit sensibly.
  expect_equal(fit$estimates$coefficients$treatment, 4, tolerance = 1)
})

test_that("bdplogit supports unordered factor covariates", {
  d <- make_factor_logit_data()
  set.seed(3007)
  fit <- bdplogit(Y ~ treatment + x + site,
    data = d$data, data0 = d$data0,
    method = "fixed", number_mcmc_beta = 3000
  )

  expect_named(
    fit$estimates$coefficients,
    c("intercept", "treatment", "x", "siteB", "siteC")
  )
  expect_true(all(is.finite(unlist(fit$estimates$coefficients))))
  expect_true(all(c("siteB", "siteC") %in% names(fit$posterior)))

  # The simulated treatment log-odds effect is 0.9.
  expect_equal(fit$estimates$coefficients$treatment, 0.9, tolerance = 0.6)
})
