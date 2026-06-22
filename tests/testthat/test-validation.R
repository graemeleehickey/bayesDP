## Tests for input validation across the main fitting functions.

test_that("bdpbinomial errors on incomplete or invalid input", {
  expect_error(bdpbinomial(N_t = 500), "Current treatment not provided")
  expect_error(
    bdpbinomial(y_t = 10, N_t = 500, y0_t = 25),
    "Historical treatment incomplete"
  )
  expect_error(
    bdpbinomial(y_t = 10, N_t = 500, discount_function = "nope"),
    "discount_function input incorrectly"
  )
})

test_that("bdpbinomial honours compare = FALSE", {
  set.seed(551)
  fit <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    method = "fixed", compare = FALSE, number_mcmc = 1e4
  )
  expect_null(fit$final)
})

test_that("bdpnormal errors on incomplete or invalid input", {
  expect_error(bdpnormal(sigma_t = 10, N_t = 50), "Current treatment not provided")
  expect_error(
    bdpnormal(mu_t = 30, sigma_t = 10, N_t = 50, mu0_t = 32),
    "Historical treatment incomplete"
  )
  expect_error(
    bdpnormal(mu_t = 30, sigma_t = 10, N_t = 50, discount_function = "nope"),
    "discount_function input incorrectly"
  )
})

test_that("bdpnormal honours compare = FALSE", {
  set.seed(552)
  fit <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 50,
    mu0_t = 32, sigma0_t = 10, N0_t = 50,
    method = "fixed", compare = FALSE, number_mcmc = 1e4
  )
  expect_null(fit$final)
})

test_that("bdpnormal and bdpbinomial run with the mc method", {
  set.seed(553)
  fit_n <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 32, sigma0_t = 10, N0_t = 250,
    method = "mc", number_mcmc = 1e4
  )
  expect_length(fit_n$posterior_treatment$alpha_discount, 1e4)
  expect_true(all(fit_n$posterior_treatment$alpha_discount >= 0 &
    fit_n$posterior_treatment$alpha_discount <= 1))

  set.seed(554)
  fit_b <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    method = "mc", number_mcmc = 1e4
  )
  expect_length(fit_b$posterior_treatment$alpha_discount, 1e4)
  expect_true(all(fit_b$posterior_treatment$alpha_discount >= 0 &
    fit_b$posterior_treatment$alpha_discount <= 1))
})

test_that("bdpbinomial and bdpnormal support the weibull discount function", {
  set.seed(555)
  fit_b <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    discount_function = "weibull", method = "fixed", number_mcmc = 1e4
  )
  expect_true(fit_b$posterior_treatment$alpha_discount >= 0 &&
    fit_b$posterior_treatment$alpha_discount <= 1)

  set.seed(556)
  fit_n <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 32, sigma0_t = 10, N0_t = 250,
    discount_function = "scaledweibull", method = "fixed", number_mcmc = 1e4
  )
  expect_true(fit_n$posterior_treatment$alpha_discount >= 0 &&
    fit_n$posterior_treatment$alpha_discount <= 1)
})

test_that("bdpsurvival errors when current data are missing", {
  expect_error(bdpsurvival(Surv(time, status) ~ 1), "Current data not input correctly")
})
