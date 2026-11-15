test_that("normal comparison agrees with the model and is invariant to outcome units", {
  for (method in c("fixed", "mc")) {
    set.seed(130)
    helper <- probability_discount(mu = 3, sigma = 2, N = 40,
      mu0 = 3.4, sigma0 = 3, N0 = 75, method = method, number_mcmc = 4000)
    set.seed(130)
    fit <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
      mu0_t = 3.4, sigma0_t = 3, N0_t = 75, method = method, number_mcmc = 4000)
    expect_equal(helper, fit$posterior_treatment$p_hat)
    set.seed(130)
    rescaled <- probability_discount(mu = 30, sigma = 20, N = 40,
      mu0 = 34, sigma0 = 30, N0 = 75, method = method, number_mcmc = 4000)
    expect_equal(helper, rescaled, tolerance = 1e-10)
  }
})

test_that("fixed binomial comparisons match integration and discount conflict in either direction", {
  for (historical_events in c(5, 20, 45)) {
    q <- integrate(function(theta) {
      dbeta(theta, 21, 81) * pbeta(theta, historical_events + 1,
        101 - historical_events, lower.tail = FALSE)
    }, 0, 1)$value
    set.seed(100 + historical_events)
    fit <- bdpbinomial(y_t = 20, N_t = 100,
      y0_t = historical_events, N0_t = 100, method = "fixed",
      alpha_max = 0.4, number_mcmc = 20000)
    expected <- 2 * min(q, 1 - q)
    expect_equal(fit$posterior_treatment$p_hat, expected, tolerance = 0.025)
    expect_equal(fit$posterior_treatment$alpha_discount, 0.4 * expected, tolerance = 0.01)
  }
})

test_that("normal example has a positive treatment-minus-control credible interval", {
  set.seed(42)
  fit <- bdpnormal(mu_t = 45, sigma_t = 10, N_t = 50,
    mu0_t = 50, sigma0_t = 10, N0_t = 50,
    mu_c = 40, sigma_c = 10, N_c = 50,
    mu0_c = 40, sigma0_c = 10, N0_c = 50, method = "fixed")
  delta <- fit$posterior_treatment$posterior_mu - fit$posterior_control$posterior_mu
  expect_gt(unname(quantile(delta, 0.025)), 0)
  expect_equal(fit$final$posterior, delta)
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "equal-tailed credible interval")
  expect_match(out, "Posterior medians")
  expect_false(grepl("alternative hypothesis", out, fixed = TRUE))
})

test_that("discount plots respect both arm limits and prescribed weights for every outcome", {
  current <- data.frame(time = c(1, 2, 3, 4, 5, 6, 1, 2, 4, 5, 6, 7),
    status = rep(c(1, 0, 1), 4), treatment = rep(c(1, 0), each = 6))
  historical <- transform(current, time = time + 0.1)
  factories <- list(
    function(...) bdpbinomial(y_t = 10, N_t = 80, y0_t = 12, N0_t = 90,
      y_c = 15, N_c = 80, y0_c = 16, N0_c = 90, ...),
    function(...) bdpnormal(mu_t = 5, sigma_t = 2, N_t = 40,
      mu0_t = 5.2, sigma0_t = 2, N0_t = 50,
      mu_c = 3, sigma_c = 2, N_c = 40, mu0_c = 3, sigma0_c = 2, N0_c = 50, ...),
    function(...) bdpsurvival(Surv(time, status) ~ treatment,
      data = current, data0 = historical, breaks = c(2.5, 4.5), ...)
  )
  for (factory in factories) {
    for (discount in c("identity", "weibull", "scaledweibull")) {
      set.seed(39)
      fit <- factory(method = "fixed", alpha_max = c(0.25, 0.65),
        weibull_scale = 0.8, weibull_shape = 3,
        discount_function = discount, number_mcmc = 500)
      p <- plot(fit, type = "discount", print = FALSE)
      endpoint <- p$data[p$data$p == 1, ]
      maxima <- c(Treatment = 0.25, Control = 0.65)[endpoint$group]
      factor <- if (discount == "weibull") 1 - exp(-(1 / 0.8)^3) else 1
      expect_equal(endpoint$alpha, unname(maxima * factor))
      expect_true(all(p$data$alpha[p$data$p == 0] == 0))
    }
    fixed <- factory(method = "fixed", alpha_max = c(0.25, 0.65),
      fix_alpha = TRUE, number_mcmc = 500)
    p <- plot(fixed, type = "discount", print = FALSE)
    expect_equal(unique(p$data$alpha[p$data$group == "Treatment"]), 0.25)
    expect_equal(unique(p$data$alpha[p$data$group == "Control"]), 0.65)
  }
})

test_that("effect plots show paired treatment-minus-control draws", {
  set.seed(75)
  fits <- list(
    bdpbinomial(y_t = 10, N_t = 50, y_c = 20, N_c = 50, number_mcmc = 2000),
    bdpnormal(mu_t = 2, sigma_t = 1, N_t = 50,
      mu_c = 10, sigma_c = 1, N_c = 50, number_mcmc = 2000)
  )
  for (fit in fits) {
    p <- plot(fit, type = "effect", print = FALSE)
    expect_s3_class(p, "ggplot")
    # Both examples have negative treatment-control effects; plotting either
    # arm's positive level instead would fail this check.
    expect_lt(p$data$effect[which.max(p$data$density)], 0)
    expect_match(p$labels$x, "treatment - control", fixed = TRUE)
    expect_equal(p$layers[[2]]$data$xintercept, 0)
  }
  one <- bdpbinomial(y_t = 10, N_t = 50, number_mcmc = 500)
  expect_error(plot(one, type = "effect", print = FALSE), "requires treatment and control")
  expect_null(plot(one, type = "discount", print = FALSE))
})

test_that("binomial summaries retain four decimals and identify posterior medians", {
  set.seed(42)
  fit <- bdpbinomial(y_t = 15, N_t = 200, y0_t = 25, N0_t = 250,
    y_c = 20, N_c = 250, y0_c = 20, N0_c = 250, method = "fixed")
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, sprintf("%.4f", median(fit$posterior_treatment$posterior)), fixed = TRUE)
  expect_match(out, "Posterior medians", fixed = TRUE)
  expect_match(out, "equal-tailed credible interval", fixed = TRUE)
})

test_that("regression covariate prior behaviour agrees with the documented method distinction", {
  set.seed(770)
  make_data <- function(n) {
    treatment <- rep(c(0, 1), length.out = n)
    x <- rnorm(n)
    data.frame(Y = 4 + treatment + 2 * x + rnorm(n), treatment = treatment, x = x)
  }
  current <- make_data(60)
  historical <- make_data(80)
  fit_model <- function(method, prior_sd, prior_mean) {
    set.seed(80)
    bdplm(Y ~ treatment + x, data = current, data0 = historical,
      method = method, prior_covariate_sd = prior_sd,
      prior_covariate_effect = prior_mean,
      number_mcmc_alpha = 1000, number_mcmc_sigmagrid = 1000,
      number_mcmc_sigma = 100, number_mcmc_beta = 2000)
  }
  mc_vague <- fit_model("mc", 1e4, 0)
  mc_nominally_strong <- fit_model("mc", 0.01, 100)
  expect_equal(mc_vague$posterior, mc_nominally_strong$posterior)
  fixed_vague <- fit_model("fixed", 1e4, 0)
  fixed_strong <- fit_model("fixed", 1e-8, 0)
  expect_gt(mean(fixed_vague$posterior$x), 1.5)
  expect_lt(abs(mean(fixed_strong$posterior$x)), 0.05)
})
