centering_data <- function(binary = FALSE) {
  current <- expand.grid(treatment = 0:1, x = c(-2, 0, 2), replicate = 1:8)
  historical <- expand.grid(treatment = 0:1, x = c(-2, 0, 2), replicate = 1:16)
  current$x <- current$x + 10 + current$treatment
  historical$x <- historical$x + 14 + 2 * historical$treatment
  current$site <- factor(ifelse(current$replicate <= 2, "A", "B"))
  historical$site <- factor(ifelse(historical$replicate <= 12, "A", "B"))
  outcome <- function(d) {
    if (binary) d$replicate %% 2 else
      2 + 0.7 * d$treatment + 0.3 * d$x + 0.6 * (d$site == "B") +
        c(-1, 0.3, -0.3, 1)[(d$replicate - 1) %% 4 + 1]
  }
  current$Y <- outcome(current)
  historical$Y <- outcome(historical)
  list(current = current, historical = historical)
}

centering_fit <- function(model, method, d, ...) {
  set.seed(9042)
  args <- list(formula = Y ~ treatment + x + site,
    data = d$current, data0 = d$historical, method = method,
    number_mcmc_alpha = 300, number_mcmc_beta = 300)
  if (model == "bdplm") args <- c(args,
    list(number_mcmc_sigmagrid = 300, number_mcmc_sigma = 50))
  do.call(get(model), c(args, list(...)))
}

test_that("centering uses the selected population for both datasets and arm priors", {
  withr::local_seed(9041)
  # Historical data have twice as many participants as current data, with
  # different continuous-covariate means and factor proportions.
  means <- list(pooled = c(x = 13.5, siteB = 5 / 12),
    current = c(x = 10.5, siteB = 3 / 4),
    historical = c(x = 15, siteB = 1 / 4))
  for (model in c("bdplm", "bdplogit")) {
    d <- centering_data(model == "bdplogit")
    for (method in if (model == "bdplm") c("fixed", "mc") else "fixed") {
      default <- centering_fit(model, method, d)
      expect_identical(default$args1$centering, "pooled")
      for (choice in names(means)) {
        fit <- centering_fit(model, method, d, centering = choice)
        expected <- means[[choice]]
        expect_identical(fit$args1$centering, choice)
        expect_equal(fit$args1$covariate_means, expected)
        if (choice == "pooled") {
          expect_identical(fit$posterior, default$posterior)
          expect_identical(fit$alpha_discount, default$alpha_discount)
        }

        raw <- rbind(d$current, d$historical)
        expect_equal(fit$args1$data$x, raw$x - expected[["x"]])
        expect_equal(fit$args1$data$siteB,
          as.numeric(raw$site == "B") - expected[["siteB"]])
        expect_equal(fit$args1$data$Y, raw$Y)
        expect_equal(fit$args1$data$treatment, raw$treatment)

        # Independently estimate the historical arm priors at the chosen
        # covariate values, and check the fit against explicitly supplied priors.
        reference_data <- transform(d$historical, control = 1 - treatment,
          x = x - expected[["x"]],
          siteB = as.numeric(site == "B") - expected[["siteB"]])
        reference_formula <- Y ~ -1 + treatment + control + x + siteB
        reference <- if (model == "bdplm") {
          lm(reference_formula, data = reference_data)
        } else {
          glm(reference_formula, data = reference_data, family = binomial())
        }
        prior <- summary(reference)$coefficients
        supplied <- centering_fit(model, method, d, centering = choice,
          prior_treatment_effect = prior["treatment", "Estimate"],
          prior_control_effect = prior["control", "Estimate"],
          prior_treatment_sd = prior["treatment", "Std. Error"],
          prior_control_sd = prior["control", "Std. Error"])
        expect_equal(fit$posterior, supplied$posterior, tolerance = 1e-7)

        # Residual reporting must use the original covariate scale for every
        # choice, including encoded factors.
        X <- model.matrix(~ treatment + x + site, d$current)
        beta <- colMeans(fit$posterior[, c("intercept", "treatment", "x", "siteB")])
        fitted <- drop(X %*% beta)
        if (model == "bdplogit") fitted <- plogis(fitted)
        residuals <- d$current$Y - fitted
        out <- withr::with_options(list(width = 200), capture.output(summary(fit)))
        observed <- scan(text = out[match("Residuals:", out) + 2L], quiet = TRUE)
        expect_equal(observed,
          unname(round(quantile(residuals, c(0, 0.25, 0.5, 0.75, 1)), 3)))

        if (choice != "pooled") {
          shifted_data <- lapply(d, function(x) transform(x, x = x + 100))
          shifted <- centering_fit(model, method, shifted_data, centering = choice)
          expect_equal(shifted$posterior$treatment, fit$posterior$treatment,
            tolerance = 1e-7)
          expect_equal(shifted$posterior$x, fit$posterior$x, tolerance = 1e-7)
          expect_equal(shifted$posterior$intercept,
            fit$posterior$intercept - 100 * fit$posterior$x, tolerance = 1e-7)
        }
      }
    }
  }
})

test_that("centering choices are independent of prescribed borrowing weights", {
  withr::local_seed(9043)
  for (model in c("bdplm", "bdplogit")) {
    d <- centering_data(model == "bdplogit")
    for (choice in c("pooled", "current", "historical")) {
      estimated <- centering_fit(model, "fixed", d, centering = choice)
      for (weights in list(0, c(0.2, 0.8))) {
        prescribed <- centering_fit(model, "fixed", d, centering = choice,
          fix_alpha = TRUE, alpha_max = weights)
        expect_identical(prescribed$args1$covariate_means, estimated$args1$covariate_means)
        expect_identical(prescribed$args1$data, estimated$args1$data)
      }
    }
  }
})

test_that("centering choices also support models without covariates", {
  withr::local_seed(9044)
  for (model in c("bdplm", "bdplogit")) {
    d <- centering_data(model == "bdplogit")
    for (method in if (model == "bdplm") c("fixed", "mc") else "fixed") {
      fits <- lapply(c("pooled", "current", "historical"), function(choice) {
        args <- list(formula = Y ~ treatment, data = d$current,
          data0 = d$historical, method = method, centering = choice,
          number_mcmc_alpha = 300, number_mcmc_beta = 300)
        if (model == "bdplm") args <- c(args,
          list(number_mcmc_sigmagrid = 300, number_mcmc_sigma = 50))
        set.seed(9045)
        do.call(get(model), args)
      })
      for (fit in fits) {
        expect_length(fit$args1$covariate_means, 0)
        expect_identical(fit$posterior, fits[[1]]$posterior)
      }
    }
  }
})

test_that("invalid centering choices are rejected before sampling", {
  withr::local_seed(9046)
  for (model in c("bdplm", "bdplogit")) {
    d <- centering_data(model == "bdplogit")
    for (choice in list("unknown", "Current", "cur", NA_character_, NULL,
        character(), c("current", "historical"), 1, TRUE)) {
      before <- .Random.seed
      expect_error(get(model)(Y ~ treatment + x, data = d$current,
        data0 = d$historical, centering = choice),
        "centering must be 'pooled', 'current', or 'historical'", fixed = TRUE)
      expect_identical(.Random.seed, before)
    }
  }
})
