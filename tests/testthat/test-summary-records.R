survival_summary_table <- function(fit) {
  withr::local_options(width = 200)
  out <- capture.output(summary(fit))
  rows <- grep("^\\s*[0-9]", out, value = TRUE)
  as.matrix(read.table(text = paste(rows, collapse = "\n"), header = FALSE))
}

regression_residual_summary <- function(fit) {
  withr::local_options(width = 200)
  out <- capture.output(summary(fit))
  scan(text = out[match("Residuals:", out) + 2L], quiet = TRUE)
}

test_that("one-arm risk sets count participants at observed times", {
  withr::local_seed(181)
  times <- c(1, 2, 2, 4, 5, 6, 7, 8)
  statuses <- list(rep(1, 8), c(1, 1, 0, 0, 1, 1, 0, 1), rep(0, 8))
  for (status in statuses) {
    current <- data.frame(time = times, status = status)
    historical <- transform(current, time = time + 0.5)
    fit <- bdpsurvival(Surv(time, status) ~ 1,
      data = current, data0 = historical, breaks = c(2, 3, 6),
      method = "fixed", number_mcmc = 1000)
    observed_times <- sort(unique(times))
    expected_counts <- cbind(observed_times,
      vapply(observed_times, function(t) sum(times >= t), integer(1)),
      vapply(observed_times, function(t) sum(times == t & status == 1), integer(1)))
    reported <- survival_summary_table(fit)
    expect_equal(unname(reported[, 1:3]), unname(expected_counts))
    expect_false(3 %in% reported[, 1]) # A cutpoint alone is not an observation.
    expected_records <- rbind(
      transform(current, treatment = 1, historical = 0),
      transform(historical, treatment = 1, historical = 1))
    expect_equal(fit$args1$participants, expected_records)

    # Survival columns still summarize the fitted hazard draws at those times.
    expected_survival <- t(vapply(observed_times, function(t) {
      draws <- as.numeric(1 - ppexp(t, fit$posterior_treatment$posterior_hazard,
        cuts = c(0, fit$args1$breaks)))
      c(median(draws), sd(draws), quantile(draws, c(0.025, 0.975)))
    }, numeric(4)))
    expect_equal(unname(reported[, 4:7]), unname(round(expected_survival, 4)))
    expect_output(summary(fit), "Posterior median survival, SD")
    out <- capture.output(print(fit))
    printed <- scan(text = grep("^\\s*[0-9]", out, value = TRUE), quiet = TRUE)
    expect_equal(printed[1:2], c(length(times), sum(status)))
  }
})

test_that("survival reports also handle renamed outcomes and no historical data", {
  withr::local_seed(182)
  current <- data.frame(followup = c(1, 2, 2, 4, 6, 8),
    event = c(1, 0, 1, 0, 1, 1))
  expect_warning(fit <- bdpsurvival(Surv(followup, event) ~ 1,
    data = current, breaks = c(2, 3, 6),
    method = "fixed", number_mcmc = 1000), "Historical treatment data missing")
  expect_equal(fit$args1$participants$time, current$followup)
  expect_equal(fit$args1$participants$status, current$event)
  expect_equal(unname(survival_summary_table(fit)[, 2]), c(6, 5, 3, 2, 1))
})

test_that("two-arm counts distinguish current and historical-only controls", {
  withr::local_seed(183)
  followup <- function(n, arm, shift = 0) {
    data.frame(time = rep(c(1, 2, 2, 4, 5, 6, 7, 8), length.out = n) + shift,
      status = rep(c(1, 0, 1, 0), length.out = n), treatment = arm)
  }
  current_t <- followup(8, 1)
  current_c <- followup(6, 0)
  historical_t <- followup(10, 1, 0.5)
  historical_c <- followup(7, 0, 0.5)
  fit <- bdpsurvival(Surv(time, status) ~ treatment,
    data = rbind(current_t, current_c), data0 = rbind(historical_t, historical_c),
    breaks = c(2, 3, 6), method = "fixed", number_mcmc = 1000)
  expect_output(summary(fit), "Current treatment: n = 8, number of events = 4")
  expect_output(summary(fit), "Current control: n = 6, number of events = 3")
  expect_output(print(fit), "Current control: n = 6, number of events = 3")
  expect_equal(nrow(fit$args1$participants), 31L)

  expect_warning(historical_control <- bdpsurvival(Surv(time, status) ~ treatment,
    data = current_t, data0 = historical_c, breaks = c(2, 3, 6),
    method = "fixed", number_mcmc = 1000), "Current control data missing")
  expect_output(summary(historical_control), "Current treatment: n = 8, number of events = 4")
  expect_output(summary(historical_control), "Historical control: n = 7, number of events = 4")
})

test_that("regression residuals agree with original-scale predictions and covariate shifts", {
  withr::local_seed(184)
  make_data <- function(n) {
    treatment <- rep(c(0, 1), length.out = n)
    x <- rnorm(n, 10, 2)
    z <- rnorm(n, -4, 1)
    site <- factor(rep(c("A", "B", "C"), length.out = n))
    data.frame(Y = 5 + 2 * treatment + x - z + as.numeric(site) + rnorm(n),
      treatment = treatment, x = x, z = z, site = site)
  }
  current <- make_data(60)
  historical <- make_data(90)
  formula <- Y ~ treatment + x + z + site
  for (method in c("fixed", "mc")) {
    fit_model <- function(current, historical) {
      set.seed(185)
      bdplm(formula, data = current, data0 = historical, method = method,
        number_mcmc_alpha = 1000, number_mcmc_sigmagrid = 500,
        number_mcmc_sigma = 100, number_mcmc_beta = 1000)
    }
    fit <- fit_model(current, historical)
    X <- model.matrix(formula, current)
    colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
    beta <- unlist(fit$estimates$coefficients[colnames(X)])
    expected_residuals <- current$Y - drop(X %*% beta)
    expected_summary <- unname(round(quantile(expected_residuals,
      c(0, 0.25, 0.5, 0.75, 1)), 3))
    expect_equal(regression_residual_summary(fit), expected_summary)

    shifted <- fit_model(transform(current, x = x + 100, z = z - 20),
      transform(historical, x = x + 100, z = z - 20))
    expect_equal(regression_residual_summary(shifted), expected_summary)
    expect_equal(shifted$posterior$treatment, fit$posterior$treatment,
      tolerance = 1e-7)

    # Reporting uses saved, named design columns, not their positions or a
    # newly evaluated formula with different factor contrasts.
    fit$args1$data <- fit$args1$data[, rev(names(fit$args1$data))]
    withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
      expect_equal(regression_residual_summary(fit), expected_summary)
    })
  }
})

test_that("logistic fits retain all coefficients and report binary response residuals", {
  withr::local_seed(187)
  make_data <- function(n) {
    treatment <- rep(c(0, 1), length.out = n)
    x <- rnorm(n, 10, 2)
    site <- factor(rep(c("A", "B", "C"), length.out = n))
    probability <- plogis(-1 + 0.6 * treatment + 0.3 * (x - 10) +
      0.2 * (site == "B") - 0.2 * (site == "C"))
    data.frame(Y = rbinom(n, 1, probability), treatment = treatment, x = x, site = site)
  }
  current <- make_data(120)
  historical <- make_data(180)
  formula <- Y ~ treatment + x + site
  fit <- bdplogit(formula, data = current, data0 = historical, method = "fixed",
    number_mcmc_alpha = 1000, number_mcmc_beta = 2000)
  X <- model.matrix(formula, current)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intercept"
  beta <- unlist(fit$estimates$coefficients[colnames(X)])
  residuals <- current$Y - plogis(drop(X %*% beta))
  expect_equal(regression_residual_summary(fit),
    unname(round(quantile(residuals, c(0, 0.25, 0.5, 0.75, 1)), 3)))
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "siteC")
  expect_false(grepl("Residual standard error", out, fixed = TRUE))
  expect_output(print(fit), "siteC")
})

test_that("older saved fits request refitting rather than reporting incorrect summaries", {
  withr::local_seed(186)
  d <- data.frame(time = 1:8, status = rep(c(1, 0), 4))
  fit <- bdpsurvival(Surv(time, status) ~ 1, data = d, data0 = d,
    breaks = c(2, 4, 6), method = "fixed", number_mcmc = 500)
  fit$args1$participants <- NULL
  expect_error(summary(fit), "participant-level records.*Refit")
  expect_error(print(fit), "participant-level records.*Refit")

  current <- data.frame(Y = rnorm(40), treatment = rep(c(0, 1), 20))
  regression <- bdplm(Y ~ treatment, data = current, data0 = current,
    method = "fixed", number_mcmc_alpha = 500, number_mcmc_sigmagrid = 500,
    number_mcmc_sigma = 100, number_mcmc_beta = 500)
  expect_length(regression$args1$covariate_means, 0)
  expect_length(regression_residual_summary(regression), 5)
  regression$args1$covariate_means <- NULL
  expect_error(summary(regression), "centering means.*Refit")
})
