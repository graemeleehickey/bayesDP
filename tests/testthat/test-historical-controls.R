test_that("historical-only binomial controls retain their Beta posterior", {
  withr::local_seed(2602)
  fit <- bdpbinomial(y_t = 8, N_t = 40, y0_c = 15, N0_c = 60,
    a0 = 2, b0 = 3, method = "fixed", number_mcmc = 20000)
  control <- fit$posterior_control
  expect_true(fit$args1$arm2)
  expect_null(control$posterior_flat)
  expect_null(control$p_hat)
  expect_null(control$alpha_discount)
  expect_identical(control$posterior, control$prior)
  expect_equal(mean(control$posterior), 17 / 65, tolerance = 0.003, scale = 1)
  expect_equal(var(control$posterior), 17 * 48 / (65^2 * 66),
    tolerance = 0.0002, scale = 1)
  expect_equal(fit$final$posterior,
    fit$posterior_treatment$posterior - control$posterior)
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "Historical control data: 15 and 60", fixed = TRUE)
  expect_false(grepl("Current control", out, fixed = TRUE))
  expect_false(grepl("Stochastic comparison", out, fixed = TRUE))
})

test_that("historical-only normal controls retain the posterior of their mean", {
  withr::local_seed(2603)
  fit <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
    mu0_c = 2, sigma0_c = 3, N0_c = 60,
    method = "fixed", number_mcmc = 20000)
  control <- fit$posterior_control
  expect_true(fit$args1$arm2)
  expect_null(control$posterior_flat_mu)
  expect_null(control$p_hat)
  expect_null(control$alpha_discount)
  expect_identical(control$posterior_mu, control$prior_mu)
  expect_equal(control$posterior_sigma2, control$prior_sigma2 / 60)
  expect_equal(mean(control$posterior_mu), 2, tolerance = 0.015, scale = 1)
  expect_equal(var(control$posterior_mu), (3^2 / 60) * 59 / 57,
    tolerance = 0.01, scale = 1)
  expect_equal(fit$final$posterior,
    fit$posterior_treatment$posterior_mu - control$posterior_mu)
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "Historical control: mu0_c = 2, sigma0_c = 3, N0_c = 60",
    fixed = TRUE)
  expect_false(grepl("Current control", out, fixed = TRUE))
  expect_false(grepl("Stochastic comparison", out, fixed = TRUE))
})
