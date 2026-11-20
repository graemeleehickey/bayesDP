survival_borrowing_data <- function(two_arm = FALSE) {
  d <- data.frame(time = c(1, 2, 2, 4, 5, 6, 7, 8),
    status = c(1, 1, 0, 0, 1, 1, 0, 1))
  if (two_arm) d <- rbind(transform(d, treatment = 1),
    transform(d, time = time * 1.1, treatment = 0))
  d
}

test_that("one-arm comparisons use the variance of the log survival difference", {
  current <- list(shape = c(2, 4, 1), rate = c(10, 12, 5),
    log_hazard = log(rbind(c(.2, .3, .1), c(.4, .1, .2), c(.1, .2, .3))))
  historical <- list(shape = c(3, 2, 5), rate = c(8, 6, 15),
    log_hazard = log(rbind(c(.3, .2, .2), c(.1, .4, .3), c(.2, .1, .4))))
  for (time in c(1, 2, 3, 5, 8)) {
    lengths <- c(min(time, 2), max(0, min(time - 2, 3)), max(0, time - 5))
    difference <- drop((exp(historical$log_hazard) - exp(current$log_hazard)) %*% lengths)
    variance <- sum(lengths^2 * (current$shape / current$rate^2 +
      historical$shape / historical$rate^2))
    expected <- 2 * pnorm(-abs(difference) / sqrt(variance))
    expect_equal(bdp_survival_comparison(current, historical, time,
      c(2, 5), FALSE, "mc"), expected, tolerance = 1e-12)
    q <- mean(difference > 0) + .5 * mean(difference == 0)
    expect_equal(bdp_survival_comparison(current, historical, time,
      c(2, 5), FALSE, "fixed"), 2 * min(q, 1 - q))
  }
})

test_that("two-arm comparisons use analytic log-Gamma inverse variances", {
  current <- list(shape = c(1, 2, 4), log_hazard = matrix(0, 5, 3))
  historical <- list(shape = c(2, 4, 1),
    log_hazard = rbind(c(0, 0, 0), c(1, 0, 0), c(1, 2, 3),
      c(-1, 2, -3), c(-2, -1, -3)))
  # At integer shapes, trigamma(k) = pi^2/6 - sum_{r=1}^{k-1} 1/r^2.
  integer_variance <- function(k) pi^2 / 6 - sum(1 / seq_len(k - 1)^2)
  variance <- vapply(current$shape, integer_variance, numeric(1)) +
    vapply(historical$shape, integer_variance, numeric(1))
  difference <- drop(historical$log_hazard %*% (1 / variance)) / sum(1 / variance)
  expected <- 2 * pnorm(-abs(difference) * sqrt(sum(1 / variance)))
  observed <- bdp_survival_comparison(current, historical, NULL, c(2, 5), TRUE, "mc")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(observed[1], 1)
  expect_gt(observed[2], observed[3])
  q <- mean(difference > 0) + .5 * mean(difference == 0)
  expect_equal(bdp_survival_comparison(current, historical, NULL,
    c(2, 5), TRUE, "fixed"), 2 * min(q, 1 - q))

  # Neither draw order, number of repeated draws, nor interval order changes weights.
  for (rows in list(c(5, 2, 1, 4, 3), rep(seq_len(5), 2))) {
    columns <- c(3, 1, 2)
    permute <- function(x) list(shape = x$shape[columns],
      log_hazard = x$log_hazard[rows, columns, drop = FALSE])
    expect_equal(bdp_survival_comparison(permute(current), permute(historical),
      NULL, c(2, 5), TRUE, "mc"), observed[rows])
  }
})

test_that("one-interval comparisons survive extreme survival probabilities", {
  current <- list(shape = 2, rate = 1, log_hazard = matrix(log(c(1, 1, 1, 1))))
  historical <- list(shape = 2, rate = 1, log_hazard = matrix(log(c(1, 2, 3, .5))))
  expected <- 2 * pnorm(-abs(c(0, 1, 2, -.5)) / 2)
  for (time in c(1e-200, 1, 1e6)) {
    expect_equal(bdp_survival_comparison(current, historical, time,
      numeric(), FALSE, "mc"), expected, tolerance = 1e-12)
    expect_equal(bdp_survival_comparison(current, historical, time,
      numeric(), FALSE, "fixed"), .75)
  }
  expect_equal(exp(-1e6 * exp(current$log_hazard)), matrix(0, 4, 1))
})

test_that("log-Gamma sampling retains small shapes without a hazard floor", {
  withr::local_seed(161)
  for (shape in c(.001, .1, 4)) {
    draws <- bdp_log_gamma(100000, shape, rate = 10)
    expect_true(all(is.finite(draws)))
    # Mean and variance of log Gamma; tolerances allow five Monte Carlo SEs.
    expect_equal(mean(draws), digamma(shape) - log(10),
      tolerance = 5 * sqrt(trigamma(shape) / length(draws)), scale = 1)
    variance_se <- sqrt((psigamma(shape, deriv = 3) + 2 * trigamma(shape)^2) / length(draws))
    expect_equal(var(draws), trigamma(shape), tolerance = 5 * variance_se, scale = 1)
    if (shape == .001) expect_true(any(exp(draws) == 0))
    if (shape == 4) {
      expect_equal(mean(exp(draws)), .4, tolerance = .005)
      expect_equal(var(exp(draws)), .04, tolerance = .002)
    }
  }
  set.seed(162)
  original <- bdp_log_gamma(1000, rep(c(.001, 2), 500), rate = 10)
  set.seed(162)
  rescaled <- bdp_log_gamma(1000, rep(c(.001, 2), 500), rate = 1e7)
  expect_equal(rescaled, original - log(1e6), tolerance = 1e-12)
})

test_that("both methods use one comparison transform and respect prescribed weights", {
  withr::local_seed(163)
  for (two_arm in c(FALSE, TRUE)) {
    d <- survival_borrowing_data(two_arm)
    formula <- if (two_arm) Surv(time, status) ~ treatment else Surv(time, status) ~ 1
    for (method in c("fixed", "mc")) {
      set.seed(163)
      fit <- bdpsurvival(formula, d, transform(d, time = time * 1.2),
        breaks = c(2, 4, 6), surv_time = 5, method = method,
        alpha_max = c(.7, .4), number_mcmc = 1000)
      set.seed(163)
      prescribed <- bdpsurvival(formula, d, transform(d, time = time * 1.2),
        breaks = c(2, 4, 6), surv_time = 5, method = method,
        alpha_max = c(.7, .4), fix_alpha = TRUE, number_mcmc = 1000)
      # Later arms use the same number of random draws even when weights differ;
      # Gamma samplers can consume different numbers of uniforms, so compare
      # p_hat exactly only for the first arm.
      expect_equal(prescribed$posterior_treatment$p_hat, fit$posterior_treatment$p_hat)
      expect_equal(prescribed$posterior_treatment$alpha_discount, .7)
      if (two_arm) expect_equal(prescribed$posterior_control$alpha_discount, .4)
      arms <- if (two_arm) c("treatment", "control") else "treatment"
      for (arm in arms) {
        draws <- fit[[paste0("posterior_", arm)]]
        maximum <- if (arm == "treatment") .7 else .4
        expect_true(all(is.finite(draws$p_hat) & draws$p_hat >= 0 & draws$p_hat <= 1))
        expect_equal(draws$alpha_discount, maximum * draws$p_hat)
        if (!two_arm) {
          # Events and exposure before each cutpoint, calculated from original rows.
          moments <- function(data) {
            starts <- c(0, 2, 4, 6)
            ends <- c(2, 4, 6, Inf)
            A <- .1 + vapply(seq_along(starts), function(j)
              sum(data$status[data$time > starts[j] & data$time <= ends[j]]), numeric(1))
            B <- .1 + vapply(seq_along(starts), function(j)
              sum(pmax(0, pmin(data$time, ends[j]) - starts[j])), numeric(1))
            list(A = A, B = B)
          }
          cur <- moments(d)
          hist <- moments(transform(d, time = time * 1.2))
          lengths <- c(2, 2, 1, 0)
          delta <- drop((draws$prior_hazard - draws$posterior_flat_hazard) %*% lengths)
          V <- sum(lengths^2 * (cur$A / cur$B^2 + hist$A / hist$B^2))
          expected <- if (method == "mc") 2 * pnorm(-abs(delta) / sqrt(V)) else
            2 * min(mean(delta > 0), mean(delta < 0))
          expect_equal(draws$p_hat, expected, tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("all discount functions are applied once to the survival comparison", {
  d <- survival_borrowing_data(TRUE)
  for (method in c("fixed", "mc")) {
    for (discount in c("identity", "weibull", "scaledweibull")) {
      withr::local_seed(164)
      fit <- bdpsurvival(Surv(time, status) ~ treatment, d, d,
        breaks = c(2, 4, 6), method = method, number_mcmc = 500,
        discount_function = discount, alpha_max = c(.8, .3),
        weibull_shape = c(2, 3), weibull_scale = c(.4, .7))
      for (arm in seq_len(2)) {
        result <- fit[[c("posterior_treatment", "posterior_control")[arm]]]
        expected <- result$p_hat
        if (discount != "identity") {
          expected <- pweibull(expected, c(2, 3)[arm], c(.4, .7)[arm])
          if (discount == "scaledweibull") expected <- expected /
            pweibull(1, c(2, 3)[arm], c(.4, .7)[arm])
        }
        expect_equal(result$alpha_discount, c(.8, .3)[arm] * expected)
      }
    }
  }
})

test_that("zero events and empty intervals retain finite comparisons and log effects", {
  withr::local_seed(165)
  for (two_arm in c(FALSE, TRUE)) {
    d <- transform(survival_borrowing_data(two_arm), status = 0)
    formula <- if (two_arm) Surv(time, status) ~ treatment else Surv(time, status) ~ 1
    for (method in c("fixed", "mc")) {
      fit <- bdpsurvival(formula, d, d, breaks = c(2, 4, 10),
        surv_time = 12, a0 = .001, method = method, number_mcmc = 1000)
      expect_true(all(is.finite(fit$posterior_treatment$p_hat)))
      expect_true(any(fit$posterior_treatment$posterior_flat_hazard == 0))
      expect_true(all(is.finite(unlist(fit$final))))
      expect_output(summary(fit), if (two_arm) "treatment" else "survival")
      if (two_arm) {
        delta <- fit$posterior_treatment$posterior_log_hazard -
          fit$posterior_control$posterior_log_hazard
        V <- apply(delta, 2, var)
        expect_equal(fit$final$posterior_loghazard, delta %*% ((1 / V) / sum(1 / V)))
      }
    }
  }
})

test_that("single intervals, boundary events, and empty final intervals are represented", {
  d <- survival_borrowing_data()
  for (breaks in list(numeric(), c(2, 8), c(2, 8, 10))) {
    withr::local_seed(166)
    fit <- bdpsurvival(Surv(time, status) ~ 1, d, d,
      breaks = breaks, surv_time = 8, number_mcmc = 500)
    expect_equal(ncol(fit$posterior_treatment$posterior_hazard), length(breaks) + 1L)
    expect_false(anyNA(fit$args1$data$interval))
    expect_equal(sum(fit$args1$S_t$status), sum(d$status))
    expect_equal(sum(fit$args1$S_t$exposure), sum(d$time))
    if (length(breaks)) expect_equal(sum(fit$args1$S_t$status[
      as.integer(fit$args1$S_t$interval) == 1]), 2)
  }
  set.seed(166)
  ordered <- bdpsurvival(Surv(time, status) ~ 1, d, d,
    breaks = c(2, 4, 6), surv_time = 5, number_mcmc = 500)
  set.seed(166)
  unordered <- bdpsurvival(Surv(time, status) ~ 1, d, d,
    breaks = c(6, 2, 4, 2), surv_time = 5, number_mcmc = 500)
  expect_equal(ordered$posterior_treatment, unordered$posterior_treatment)

  d2 <- survival_borrowing_data(TRUE)
  for (method in c("fixed", "mc")) {
    fit <- bdpsurvival(Surv(time, status) ~ treatment, d2,
      transform(d2, time = time * 1.2), breaks = numeric(),
      method = method, number_mcmc = 500)
    expect_equal(ncol(fit$posterior_treatment$posterior_hazard), 1L)
    expect_true(all(as.integer(fit$args1$data$interval) == 1L))
    expect_equal(nrow(fit$args1$data), 2L * nrow(d2))
    expect_true(all(is.finite(fit$final$posterior_loghazard)))
  }
})

test_that("survival borrowing and treatment effects are invariant to time units", {
  for (two_arm in c(FALSE, TRUE)) {
    d <- survival_borrowing_data(two_arm)
    formula <- if (two_arm) Surv(time, status) ~ treatment else Surv(time, status) ~ 1
    for (method in c("fixed", "mc")) {
      fit_units <- function(scale) {
        set.seed(167)
        bdpsurvival(formula, transform(d, time = time * scale),
          transform(d, time = time * 1.2 * scale), breaks = c(2, 4, 6) * scale,
          surv_time = 5 * scale, b0 = .1 * scale, method = method, number_mcmc = 1000)
      }
      original <- fit_units(1)
      for (scale in c(1e-6, 365.25, 1e6)) {
        rescaled <- fit_units(scale)
        expect_equal(rescaled$posterior_treatment$p_hat, original$posterior_treatment$p_hat,
          tolerance = 1e-10)
        expect_equal(rescaled$posterior_treatment$posterior_log_hazard,
          original$posterior_treatment$posterior_log_hazard - log(scale), tolerance = 1e-10)
        expect_equal(rescaled$final, original$final, tolerance = 1e-10)
        if (two_arm) expect_equal(rescaled$posterior_control$p_hat,
          original$posterior_control$p_hat, tolerance = 1e-10)
      }
    }
  }
})

test_that("historical-posterior augmentation remains unchanged pending issue 17", {
  withr::local_seed(168)
  d <- survival_borrowing_data()
  d0 <- transform(d, time = time * 1.2)
  fit <- bdpsurvival(Surv(time, status) ~ 1, d, d0, breaks = numeric(),
    surv_time = 5, a0 = 2, b0 = 3, fix_alpha = TRUE, alpha_max = .4,
    method = "fixed", number_mcmc = 50000)
  A <- 2 + sum(d$status) + .4 * (2 + sum(d0$status))
  B <- 3 + sum(d$time) + .4 * (3 + sum(d0$time))
  draws <- as.numeric(fit$posterior_treatment$posterior_hazard)
  expect_equal(mean(draws), A / B, tolerance = .002)
  expect_equal(var(draws), A / B^2, tolerance = .0001)
})

test_that("invalid survival comparison settings give informative errors", {
  d <- survival_borrowing_data()
  fit <- function(...) bdpsurvival(Surv(time, status) ~ 1, d, d,
    breaks = numeric(), ...)
  for (value in list(0, -1, Inf, NA_real_, c(.1, .2))) {
    expect_error(fit(a0 = value), "a0 must")
    expect_error(fit(b0 = value), "b0 must")
    expect_error(fit(surv_time = value), "surv_time must")
  }
  for (value in list(1, 2.5, Inf, NA_real_)) expect_error(fit(number_mcmc = value), "number_mcmc must")
  expect_error(fit(method = "other"), "Unrecognized method")
  expect_error(fit(alpha_max = 1.1), "alpha_max must")
  expect_error(fit(discount_function = "weibull", weibull_shape = 0), "weibull_shape must")
  expect_error(fit(discount_function = "weibull", weibull_scale = Inf), "weibull_scale must")
})
