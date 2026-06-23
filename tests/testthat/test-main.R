## Tests for bayesDP fitting functions.
##
## Where possible, computed posterior values are pinned against closed-form
## (conjugate) results. For the augmented one-arm fits we use fix_alpha = TRUE
## and method = "fixed" so the historical-data weight (alpha) is known exactly,
## which makes the augmented posterior analytically tractable.
##
## All plot() calls pass an explicit `type` so the tests never trigger the
## interactive "Hit <Return> to see next plot" prompt.

## A large MCMC sample is used so Monte Carlo error is small enough to pin
## posterior means to a few decimal places.
n_mcmc <- 1e5

################################################################################
# bdpnormal
################################################################################

test_that("bdpnormal one-arm augmented posterior matches conjugate mean", {
  set.seed(7724)
  fit <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 50, sigma0_t = 5, N0_t = 250,
    fix_alpha = TRUE, method = "fixed", number_mcmc = n_mcmc
  )

  # With fix_alpha = TRUE the historical weight is exactly alpha_max = 1.
  expect_equal(fit$posterior_treatment$alpha_discount, 1)

  # Precision-weighted (conjugate) posterior mean of mu, with effective
  # historical sample size N0 * alpha = 250.
  prec_curr <- 250 / 10^2
  prec_hist <- 250 / 5^2
  analytic_mu <- (prec_curr * 30 + prec_hist * 50) / (prec_curr + prec_hist)

  expect_equal(mean(fit$posterior_treatment$posterior_mu), analytic_mu,
    tolerance = 0.02
  )

  # Structural invariants
  expect_false(fit$args1$arm2)
  expect_true(fit$args1$compare)
  expect_length(fit$posterior_treatment$posterior_mu, n_mcmc)
  expect_length(fit$final$posterior, n_mcmc)
})

test_that("bdpnormal flat-prior draw is the conjugate posterior, not the predictive", {
  # With only current data and a non-informative prior, the marginal posterior
  # of mu (posterior_flat_mu) is Student-t_{N-1}(mu_t, sigma_t^2 / N_t). Its
  # variance is (sigma^2 / N) * (N - 1) / (N - 3). This pins the construction
  # against the conjugate posterior and rules out a posterior-predictive draw
  # (which would instead scale by sigma^2 * (1 + 1/N), ~50x larger here).
  mu_t <- 30
  sigma_t <- 10
  N_t <- 50

  set.seed(6691)
  fit <- bdpnormal(
    mu_t = mu_t, sigma_t = sigma_t, N_t = N_t,
    method = "fixed", number_mcmc = 5e5
  )

  post_flat_mu <- fit$posterior_treatment$posterior_flat_mu

  conjugate_var <- (sigma_t^2 / N_t) * (N_t - 1) / (N_t - 3)
  predictive_var <- (sigma_t^2 * (1 + 1 / N_t)) * (N_t - 1) / (N_t - 3)

  expect_equal(mean(post_flat_mu), mu_t, tolerance = 0.05)
  expect_equal(var(post_flat_mu), conjugate_var, tolerance = 0.05 * conjugate_var)
  # Sanity check that we are nowhere near the posterior-predictive variance.
  expect_lt(var(post_flat_mu), predictive_var / 10)
})

test_that("bdpnormal two-arm summary, plot, and comparison object are correct", {
  set.seed(2218)
  fit <- bdpnormal(
    mu_t = 30, sigma_t = 10, N_t = 250,
    mu0_t = 50, sigma0_t = 5, N0_t = 250,
    mu_c = 25, sigma_c = 10, N_c = 250,
    mu0_c = 50, sigma0_c = 5, N0_c = 250,
    method = "fixed", number_mcmc = n_mcmc
  )

  expect_true(fit$args1$arm2)

  # The two-arm comparison object is treatment - control.
  expect_equal(
    fit$final$posterior,
    fit$posterior_treatment$posterior_mu - fit$posterior_control$posterior_mu
  )

  # summary() should report a two-arm analysis (guards the args1$arm2 bug).
  expect_output(summary(fit), "Two-armed bdp normal")

  # plot() with an explicit type returns a ggplot without prompting.
  expect_s3_class(plot(fit, type = "density", print = FALSE), "ggplot")
})

################################################################################
# bdpbinomial
################################################################################

test_that("bdpbinomial one-arm augmented posterior matches conjugate mean", {
  set.seed(513)
  fit <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    fix_alpha = TRUE, method = "fixed", number_mcmc = n_mcmc
  )

  expect_equal(fit$posterior_treatment$alpha_discount, 1)

  # Augmented posterior: Beta(y + (y0/N0)*N0*alpha + a0,
  #                           (N - y) + (N0 - (y0/N0)*N0)*alpha + b0)
  # with a0 = b0 = 1, alpha = 1.
  a_post <- 10 + (25 / 250) * 250 + 1
  b_post <- (500 - 10) + (250 - (25 / 250) * 250) + 1
  analytic_mean <- a_post / (a_post + b_post)

  expect_equal(mean(fit$posterior_treatment$posterior), analytic_mean,
    tolerance = 1e-3
  )

  expect_false(fit$args1$arm2)
  expect_true(fit$args1$compare)
  expect_true(all(fit$posterior_treatment$posterior > 0 &
    fit$posterior_treatment$posterior < 1))
})

test_that("bdpbinomial two-arm comparison, summary, print, and plot work", {
  set.seed(4406)
  fit <- bdpbinomial(
    y_t = 10, N_t = 500, y0_t = 25, N0_t = 250,
    y_c = 8, N_c = 500, y0_c = 20, N0_c = 250,
    method = "fixed", number_mcmc = n_mcmc
  )

  expect_true(fit$args1$arm2)
  expect_equal(
    fit$final$posterior,
    fit$posterior_treatment$posterior - fit$posterior_control$posterior
  )

  # Guards the args1$arm2 summary bug.
  expect_output(summary(fit), "Two-armed bdp binomial")
  expect_output(print(fit), "Two-armed bdp binomial")
  expect_s3_class(plot(fit, type = "density", print = FALSE), "ggplot")
})

################################################################################
# bdpsurvival
################################################################################

test_that("bdpsurvival one-arm survival probabilities are valid", {
  set.seed(9931)
  df_ <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 20))
  df_$status <- ifelse(df_$time < df_$status, 1, 0)
  df0 <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 10))
  df0$status <- ifelse(df0$time < df0$status, 1, 0)

  fit <- bdpsurvival(Surv(time, status) ~ 1,
    data = df_, data0 = df0,
    surv_time = 5, fix_alpha = TRUE, method = "fixed",
    number_mcmc = n_mcmc
  )

  # fix_alpha = TRUE pins the weight at alpha_max = 1.
  expect_equal(fit$posterior_treatment$alpha_discount, 1)

  ps <- fit$posterior_treatment$posterior_survival
  expect_length(ps, n_mcmc)
  expect_true(all(ps >= 0 & ps <= 1))

  expect_false(fit$args1$arm2)
  expect_s3_class(plot(fit, type = "survival", print = FALSE), "ggplot")
})

test_that("bdpsurvival default surv_time uses combined current + historical times", {
  set.seed(3185)
  df_ <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 20))
  df_$status <- ifelse(df_$time < df_$status, 1, 0)
  df0 <- data.frame(status = rexp(50, rate = 1 / 30), time = rexp(50, rate = 1 / 10))
  df0$status <- ifelse(df0$time < df0$status, 1, 0)

  # surv_time left NULL: should be the median of pooled current + historical
  # times (guards the invalid Y[, 0] index bug).
  fit <- bdpsurvival(Surv(time, status) ~ 1,
    data = df_, data0 = df0,
    fix_alpha = TRUE, method = "fixed", number_mcmc = 1e4
  )

  expect_equal(fit$args1$surv_time, median(c(df_$time, df0$time)))
})

test_that("bdpsurvival two-arm hazard-ratio summary runs and is finite", {
  set.seed(6017)
  time <- c(
    rexp(50, rate = 1 / 20), # Current treatment
    rexp(50, rate = 1 / 10), # Current control
    rexp(50, rate = 1 / 30), # Historical treatment
    rexp(50, rate = 1 / 5) # Historical control
  )
  status <- rexp(200, rate = 1 / 40)
  status <- ifelse(time < status, 1, 0)

  d2 <- data.frame(
    status = status, time = time,
    historical = c(rep(0, 100), rep(1, 100)),
    treatment = c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 50))
  )

  # The historical-data warning here is expected given this data layout and is
  # orthogonal to the behaviour under test.
  fit <- suppressWarnings(
    bdpsurvival(Surv(time, status) ~ historical + treatment,
      data = d2, method = "fixed", number_mcmc = n_mcmc
    )
  )

  expect_true(fit$args1$arm2)
  expect_true(all(is.finite(fit$final$posterior_loghazard)))

  # Two-arm summary must not error (guards the undefined data_c bug).
  expect_output(summary(fit), "Two-armed bdp survival")
})

################################################################################
# bdplm
################################################################################

test_that("bdplm recovers the treatment effect and stores discount weights", {
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

  # Treatment effect should be in a plausible neighbourhood of the truth (31).
  expect_equal(fit$estimates$coefficients$treatment, 31, tolerance = 3)

  # Discount weights are valid probabilities in [0, 1].
  alpha <- unlist(fit$alpha_discount)
  expect_true(all(alpha >= 0 & alpha <= 1))

  expect_output(summary(fit), "Coefficients:")
})
