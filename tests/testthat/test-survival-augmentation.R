survival_augmentation_data <- function() {
  list(
    current = data.frame(time = c(1, 3, 4, 5), status = c(1, 1, 0, 1)),
    historical = data.frame(time = c(1.5, 2, 3, 6), status = c(0, 1, 0, 1)),
    control = data.frame(time = c(.5, 1, 2, 4), status = c(1, 0, 1, 0)),
    control0 = data.frame(time = c(1, 2, 5, 7), status = c(0, 1, 1, 0)))
}

# Check independent Gamma draws, allowing each draw its own shape and rate.
# Conditional on the weights, Monte Carlo SEs follow from Gamma moments.
expect_survival_gamma <- function(draws, shape, rate) {
  draws <- as.numeric(draws)
  n <- length(draws)
  shape <- rep_len(shape, n)
  rate <- rep_len(rate, n)
  means <- shape / rate
  variances <- shape / rate^2
  second_moments <- shape * (shape + 1) / rate^2
  squared_draw_variances <- 2 * shape * (shape + 1) * (2 * shape + 3) / rate^4
  mean_se <- sqrt(mean(variances) / n)
  second_moment_se <- sqrt(mean(squared_draw_variances) / n)
  expect_true(all(is.finite(draws)))
  expect_equal(mean(draws), mean(means), tolerance = 6 * mean_se, scale = 1)
  expect_equal(mean(draws^2), mean(second_moments),
    tolerance = 6 * second_moment_se, scale = 1)
  # For a mixture, variance includes variation between conditional means.
  mixture_variance <- mean(variances) + mean((means - mean(means))^2)
  expect_equal(mean((draws - mean(draws))^2), mixture_variance,
    tolerance = 6 * second_moment_se + 12 * mean(means) * mean_se + 36 * mean_se^2,
    scale = 1)
  # The conditional Gamma CDF also detects incorrect pairing of vector weights.
  probabilities <- pgamma(draws, shape = shape, rate = rate)
  expect_equal(mean(probabilities), .5, tolerance = 6 / sqrt(12 * n), scale = 1)
}

test_that("prescribed survival weights discount event counts and exposure only", {
  withr::local_seed(171)
  d <- survival_augmentation_data()
  # At cuts 2 and 8, the original records give these interval totals.
  D <- c(1, 2, 0)
  T <- c(7, 6, 0)
  D0 <- c(1, 1, 0)
  T0 <- c(7.5, 5, 0)
  for (alpha in c(0, .4, 1)) {
    fit <- bdpsurvival(Surv(time, status) ~ 1, d$current, d$historical,
      breaks = c(2, 8), surv_time = 5, a0 = 2, b0 = 3,
      fix_alpha = TRUE, alpha_max = alpha, method = "fixed", number_mcmc = 20000)
    shape <- 2 + D + alpha * D0
    rate <- 3 + T + alpha * T0
    expect_equal(fit$posterior_treatment$alpha_discount, alpha)
    for (j in seq_along(shape)) {
      expect_survival_gamma(fit$posterior_treatment$posterior_hazard[, j], shape[j], rate[j])
    }
    # The Laplace transform gives the exact mean survival under independent hazards.
    lengths <- c(2, 3, 0)
    expected_survival <- prod((rate / (rate + lengths))^shape)
    expected_squared_survival <- prod((rate / (rate + 2 * lengths))^shape)
    expect_equal(mean(fit$final$posterior_survival), expected_survival,
      tolerance = 6 * sqrt((expected_squared_survival - expected_survival^2) / 20000),
      scale = 1)

    if (alpha %in% c(0, 1)) {
      reference_data <- if (alpha == 0) d$current else rbind(d$current, d$historical)
      expect_warning(reference <- bdpsurvival(Surv(time, status) ~ 1, reference_data,
        breaks = c(2, 8), surv_time = 5, a0 = 2, b0 = 3,
        method = "fixed", number_mcmc = 20000), "Historical treatment data missing")
      # Compare distributions, not identical draws from different RNG sequences.
      for (j in seq_along(shape)) {
        expect_survival_gamma(reference$posterior_treatment$posterior_hazard[, j],
          shape[j], rate[j])
        difference <- mean(fit$posterior_treatment$posterior_hazard[, j]) -
          mean(reference$posterior_treatment$posterior_hazard[, j])
        expect_equal(difference, 0,
          tolerance = 6 * sqrt(2 * shape[j] / rate[j]^2 / 20000), scale = 1)
      }
    }
  }
})

test_that("two-arm augmentation respects separate prescribed weights under both methods", {
  withr::local_seed(172)
  d <- survival_augmentation_data()
  current <- rbind(transform(d$current, treatment = 1), transform(d$control, treatment = 0))
  historical <- rbind(transform(d$historical, treatment = 1), transform(d$control0, treatment = 0))
  D <- list(c(1, 2, 0), c(2, 0, 0))
  T <- list(c(7, 6, 0), c(5.5, 2, 0))
  D0 <- list(c(1, 1, 0), c(1, 1, 0))
  T0 <- list(c(7.5, 5, 0), c(7, 8, 0))
  alpha <- c(.25, .8)
  for (method in c("fixed", "mc")) {
    fit <- bdpsurvival(Surv(time, status) ~ treatment, current, historical,
      breaks = c(2, 8), a0 = 2, b0 = 3, fix_alpha = TRUE,
      alpha_max = alpha, method = method, number_mcmc = 20000)
    for (arm in seq_len(2)) {
      result <- fit[[c("posterior_treatment", "posterior_control")[arm]]]
      expect_equal(result$alpha_discount, alpha[arm])
      for (j in seq_len(3)) {
        expect_survival_gamma(result$posterior_hazard[, j],
          2 + D[[arm]][j] + alpha[arm] * D0[[arm]][j],
          3 + T[[arm]][j] + alpha[arm] * T0[[arm]][j])
      }
    }
  }
})

test_that("estimated survival weights give the corresponding Gamma or mixture moments", {
  withr::local_seed(173)
  d <- survival_augmentation_data()
  D <- list(c(1, 2, 0), c(2, 0, 0))
  T <- list(c(7, 6, 0), c(5.5, 2, 0))
  D0 <- list(c(1, 1, 0), c(1, 1, 0))
  T0 <- list(c(7.5, 5, 0), c(7, 8, 0))
  for (two_arm in c(FALSE, TRUE)) {
    formula <- if (two_arm) Surv(time, status) ~ treatment else Surv(time, status) ~ 1
    current <- if (two_arm) rbind(transform(d$current, treatment = 1),
      transform(d$control, treatment = 0)) else d$current
    historical <- if (two_arm) rbind(transform(d$historical, treatment = 1),
      transform(d$control0, treatment = 0)) else d$historical
    for (method in c("fixed", "mc")) {
      fit <- bdpsurvival(formula, current, historical, breaks = c(2, 8),
        surv_time = 5, a0 = 2, b0 = 3, alpha_max = c(.7, .9),
        fix_alpha = FALSE, method = method, number_mcmc = 20000)
      for (arm in seq_len(if (two_arm) 2 else 1)) {
        result <- fit[[c("posterior_treatment", "posterior_control")[arm]]]
        alpha <- result$alpha_discount
        if (method == "mc") {
          expect_length(alpha, 20000)
          expect_gt(sd(alpha), .05)
        } else {
          expect_length(alpha, 1)
        }
        for (j in seq_len(3)) {
          expect_survival_gamma(result$posterior_hazard[, j],
            2 + D[[arm]][j] + alpha * D0[[arm]][j],
            3 + T[[arm]][j] + alpha * T0[[arm]][j])
        }
      }
    }
  }
})

test_that("zero-event and unobserved intervals include the baseline prior once", {
  withr::local_seed(174)
  d <- survival_augmentation_data()
  fit <- bdpsurvival(Surv(time, status) ~ 1,
    transform(d$current, status = 0), transform(d$historical, status = 0),
    breaks = c(2, 8), surv_time = 10, fix_alpha = TRUE,
    alpha_max = 1, method = "mc", number_mcmc = 20000)
  # Default Gamma(.1, .1); the empty final interval must retain that prior.
  for (j in seq_len(3)) {
    expect_survival_gamma(fit$posterior_treatment$posterior_hazard[, j],
      .1, .1 + c(7, 6, 0)[j] + c(7.5, 5, 0)[j])
  }
  expect_true(all(is.finite(fit$final$posterior_survival)))
})

test_that("source totals remain distinct from the baseline prior", {
  withr::local_seed(175)
  records <- data.frame(interval = factor(c(1, 1, 1, 1, 2, 2, 2), levels = 1:3),
    status = c(1, 0, 0, 0, 1, 0, 1), exposure = c(1, 2, 2, 2, 1, 2, 3))
  for (prior in c(.1, 1e20)) {
    source <- bdp_survival_source(records, nInt = 3, a0 = prior, b0 = prior, n = 2)
    expect_equal(source$events, c(1, 2, 0))
    expect_equal(source$exposure, c(7, 6, 0))
    expect_equal(source$shape, prior + c(1, 2, 0))
    expect_equal(source$rate, prior + c(7, 6, 0))
  }
})
