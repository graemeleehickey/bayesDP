test_that("summary-data models identify every missing component of an arm", {
  models <- list(
    list(fit = bdpbinomial, current = list(y_t = 8, N_t = 40),
      arms = list(
        list(y_t = 8, N_t = 40), list(y0_t = 10, N0_t = 50),
        list(y_c = 12, N_c = 40), list(y0_c = 15, N0_c = 50)
      )),
    list(fit = bdpnormal, current = list(mu_t = 3, sigma_t = 2, N_t = 40),
      arms = list(
        list(mu_t = 3, sigma_t = 2, N_t = 40),
        list(mu0_t = 4, sigma0_t = 2, N0_t = 50),
        list(mu_c = 2, sigma_c = 2, N_c = 40),
        list(mu0_c = 2, sigma0_c = 2, N0_c = 50)
      ))
  )
  messages <- c("Current treatment", "Historical treatment",
    "Current control", "Historical Control")
  for (model in models) {
    for (arm in seq_along(model$arms)) {
      for (missing in names(model$arms[[arm]])) {
        incomplete <- model$arms[[arm]]
        incomplete[[missing]] <- NULL
        args <- if (arm == 1L) incomplete else c(model$current, incomplete)
        expect_output(
          expect_error(do.call(model$fit, args), messages[arm]),
          paste0(missing, " missing"), fixed = TRUE
        )
      }
    }
  }
})

test_that("regression models reject incompatible formulas and treatment coding", {
  current <- data.frame(Y = rep(0:1, 6), treatment = rep(0:1, each = 6),
    x = seq_len(12), site = factor(rep(c("A", "B"), 6)))
  historical <- current
  for (fit in list(bdplm, bdplogit)) {
    expect_error(fit(Y ~ treatment, data0 = historical), "Current data")
    expect_error(fit(Y ~ treatment, data = current), "Historical data")
    for (formula in list(Y ~ x, Y ~ treatment - 1, Y ~ 0)) {
      expect_error(fit(formula, data = current, data0 = historical,
        method = "fixed"), "Intercept and treatment must be present")
    }
    expect_error(fit(Y ~ treatment, data = transform(current, treatment = treatment + 1),
      data0 = historical, method = "fixed"), "Treatment input has wrong levels")
    expect_error(fit(Y ~ treatment, data = current,
      data0 = transform(historical, treatment = treatment + 1),
      method = "fixed"), "Treatment input has wrong levels")
    historical$site <- factor(rep(c("A", "C"), 6))
    expect_error(fit(Y ~ treatment + site, data = current, data0 = historical,
      method = "fixed"), "Column names in the historical data must match")
    expect_error(fit(Y ~ treatment, data = current, data0 = historical,
      method = "fixed", discount_function = "unknown"),
      "discount_function input incorrectly")
  }
})

test_that("linear regression rejects incompatible Monte Carlo settings", {
  d <- data.frame(Y = seq_len(12), treatment = rep(0:1, each = 6))
  expect_error(bdplm(Y ~ treatment, data = d, data0 = d, method = "mc",
    number_mcmc_alpha = 100, number_mcmc_sigmagrid = 101),
    "number_mcmc_sigmagrid must be <= number_mcmc_alpha", fixed = TRUE)
  expect_error(bdplm(Y ~ treatment, data = d, data0 = d, method = "mc",
    fix_alpha = TRUE), "mc method not possible with fix_alpha")
})

test_that("survival models reject invalid cutpoints and borrowing switches", {
  withr::local_seed(2601)
  d <- data.frame(time = seq_len(6), status = c(1, 0, 1, 1, 0, 1))
  for (breaks in list(-1, NA_real_, Inf, "2")) {
    expect_error(bdpsurvival(Surv(time, status) ~ 1, data = d, data0 = d,
      breaks = breaks), "breaks must contain finite, non-negative cutpoints")
  }
  for (fixed in list(NA, 1, c(TRUE, FALSE))) {
    expect_error(bdpsurvival(Surv(time, status) ~ 1, data = d, data0 = d,
      fix_alpha = fixed), "fix_alpha must be TRUE or FALSE")
  }
  expect_error(bdpsurvival(time ~ 1, data = d), "Data input incorrectly")
  expect_error(bdpsurvival(Surv(time, status) ~ 1, data = d, data0 = d,
    discount_function = "unknown"), "discount_function input incorrectly")

  # A zero cutpoint is redundant; removing it must preserve the fitted draws.
  expect_warning(with_zero <- bdpsurvival(Surv(time, status) ~ 1,
    data = d, data0 = d, breaks = c(0, 2, 4), number_mcmc = 200),
    "zero value was removed")
  set.seed(2601)
  positive <- bdpsurvival(Surv(time, status) ~ 1,
    data = d, data0 = d, breaks = c(2, 4), number_mcmc = 200)
  expect_equal(with_zero$posterior_treatment, positive$posterior_treatment)
  expect_equal(with_zero$args1$breaks, c(2, 4))
})
