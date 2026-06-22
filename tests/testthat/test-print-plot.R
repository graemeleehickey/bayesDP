## Tests for print methods and the remaining plot branches.

test_that("print methods for bdpnormal and bdpbinomial delegate to summary", {
  set.seed(331)
  fit_n <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 50, sigma0_t = 5, N0_t = 250,
    method = "fixed", number_mcmc = 1e4
  )
  expect_output(print(fit_n), "bdp normal")

  set.seed(332)
  fit_b <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    method = "fixed", number_mcmc = 1e4
  )
  expect_output(print(fit_b), "bdp binomial")
})

test_that("one-arm bdpsurvival print produces a survival-probability table", {
  set.seed(9931)
  df_ <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 20))
  df_$status <- ifelse(df_$time < df_$status, 1, 0)
  df0 <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 10))
  df0$status <- ifelse(df0$time < df0$status, 1, 0)

  fit <- bdpsurvival(Surv(time, status) ~ 1,
    data = df_, data0 = df0,
    surv_time = 5, method = "fixed", number_mcmc = 1e4
  )

  out <- capture.output(print(fit))
  expect_true(any(grepl("One-armed bdp survival", out)))
  expect_true(any(grepl("surv_time", out)))
})

test_that("bdpnormal plot returns each plot type as a ggplot object", {
  set.seed(771)
  fit <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 50, sigma0_t = 5, N0_t = 250,
    method = "fixed", number_mcmc = 1e4
  )

  expect_s3_class(plot(fit, type = "posteriors", print = FALSE), "ggplot")
  expect_s3_class(plot(fit, type = "density", print = FALSE), "ggplot")
  expect_s3_class(plot(fit, type = "discount", print = FALSE), "ggplot")
})

test_that("bdpbinomial plot discount branch returns a ggplot object", {
  set.seed(772)
  fit <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    method = "fixed", number_mcmc = 1e4
  )
  expect_s3_class(plot(fit, type = "discount", print = FALSE), "ggplot")
})

test_that("bdpsurvival plot discount branch returns a ggplot object", {
  set.seed(773)
  df_ <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 20))
  df_$status <- ifelse(df_$time < df_$status, 1, 0)
  df0 <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 10))
  df0$status <- ifelse(df0$time < df0$status, 1, 0)

  fit <- bdpsurvival(Surv(time, status) ~ 1,
    data = df_, data0 = df0,
    surv_time = 5, method = "fixed", number_mcmc = 1e4
  )
  expect_s3_class(plot(fit, type = "discount", print = FALSE), "ggplot")
})

test_that("bdplm print and summary produce coefficient output", {
  set.seed(8852)
  n_t <- 30
  n_c <- 30
  n_t0 <- 80
  n_c0 <- 80
  treatment <- c(rep(1, n_t), rep(0, n_c))
  treatment0 <- c(rep(1, n_t0), rep(0, n_c0))
  x <- rnorm(n_t + n_c, 1, 5)
  x0 <- rnorm(n_t0 + n_c0, 1, 5)
  Y <- 10 + 31 * treatment + x * 3 + rnorm(n_t + n_c, 0, 5)
  Y0 <- 10 + 30 * treatment0 + x0 * 3 + rnorm(n_t0 + n_c0, 0, 5)
  df_ <- data.frame(Y = Y, treatment = treatment, x = x)
  df0 <- data.frame(Y = Y0, treatment = treatment0, x = x0)

  fit <- bdplm(Y ~ treatment + x, data = df_, data0 = df0, method = "fixed")

  expect_output(print(fit), "Coefficients:")
  expect_output(print(fit), "Discount function value")
})

test_that("two-arm bdpsurvival print delegates to summary", {
  set.seed(6017)
  time <- c(
    rexp(50, rate = 1 / 20), rexp(50, rate = 1 / 10),
    rexp(50, rate = 1 / 30), rexp(50, rate = 1 / 5)
  )
  status <- rexp(200, rate = 1 / 40)
  status <- ifelse(time < status, 1, 0)
  d2 <- data.frame(
    status = status, time = time,
    historical = c(rep(0, 100), rep(1, 100)),
    treatment = c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 50))
  )

  fit <- suppressWarnings(
    bdpsurvival(Surv(time, status) ~ historical + treatment,
      data = d2, method = "fixed", number_mcmc = 1e4
    )
  )
  expect_output(print(fit), "Two-armed bdp survival")
})
