extractor_fits <- make_tidier_fits()

test_that("regression extractors preserve joint posterior uncertainty and coefficient order", {
  for (fit in extractor_fits[c("lm_fixed", "lm_mc", "logit")]) {
    saved <- fit
    draws <- fit$posterior[, setdiff(names(fit$posterior), "sigma"), drop = FALSE]
    names(draws)[names(draws) == "intercept"] <- "(Intercept)"
    seed <- .Random.seed
    expect_equal(coef(fit), colMeans(draws))
    expect_identical(stats::coefficients(fit), coef(fit))
    expect_equal(vcov(fit), cov(draws))
    expect_equal(sqrt(diag(vcov(fit))), setNames(tidy(fit)$std.error, names(draws)))
    expect_equal(vcov(fit, complete = FALSE), vcov(fit))
    expect_equal(confint(fit), t(vapply(draws, quantile, numeric(2), probs = c(0.025, 0.975))))
    for (parm in list("treatment", c("siteB", "(Intercept)"), c(4L, 1L))) {
      expected <- t(vapply(draws[, parm, drop = FALSE], quantile, numeric(2), probs = c(0.1, 0.9)))
      expect_equal(confint(fit, parm, level = 0.8), expected)
    }
    expect_identical(.Random.seed, seed)
    expect_identical(fit, saved)
  }
  fit <- extractor_fits$lm_fixed
  fit$args1$covariate_means <- NULL
  expect_equal(coef(fit), coef(extractor_fits$lm_fixed))
  expect_equal(vcov(fit), vcov(extractor_fits$lm_fixed))
})

test_that("extractor arguments fail explicitly when unsupported or invalid", {
  fit <- extractor_fits$lm_fixed
  for (parm in list("sigma", "unknown", c("treatment", NA), numeric(), 0, -1, 1.5, Inf, NA_real_, 99, TRUE)) {
    expect_error(confint(fit, parm), "parm must")
  }
  for (level in list(NULL, NA_real_, Inf, 0, 1, c(0.8, 0.9), "0.9")) {
    expect_error(confint(fit, level = level), "level must")
  }
  expect_error(vcov(fit, complete = NA), "complete must")
  expect_error(coef(fit, exponentiate = TRUE), "Unsupported argument")
  expect_error(confint(fit, conf.level = 0.8), "Unsupported argument")
  expect_error(nobs(fit, historical = TRUE), "Unsupported argument")
})

test_that("nobs counts current participants rather than historical or split records", {
  withr::local_seed(3050)
  binomial <- bdpbinomial(y_t = 8, N_t = 40, y0_c = 12, N0_c = 60,
    method = "fixed", number_mcmc = 100)
  normal <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
    mu0_t = 2.8, sigma0_t = 2, N0_t = 60,
    method = "fixed", number_mcmc = 100)
  expect_equal(nobs(binomial), 40)
  expect_equal(nobs(normal), 40)
  for (fit in extractor_fits) expect_identical(nobs(fit), glance(fit)$nobs)
  expect_equal(nobs(extractor_fits$lm_mc), 48)
  expect_equal(nobs(extractor_fits$logit), 48)
  expect_equal(nobs(extractor_fits$survival_one), 6)
  expect_equal(nobs(extractor_fits$survival_two), 12)
  fit <- extractor_fits$survival_two
  fit$args1$participants <- NULL
  expect_identical(nobs(fit), NA_real_)
  fit <- extractor_fits$lm_fixed
  fit$args1$data <- NULL
  expect_identical(nobs(fit), NA_real_)
})

test_that("posterior conversion preserves arm draws and paired contrasts in all scenarios", {
  skip_if_not_installed("posterior")
  withr::local_seed(3051)
  scenarios <- list(current_only = list(N_t = 40),
    one_arm = list(N_t = 40, N0_t = 60),
    two_arm = list(N_t = 40, N0_t = 60, N_c = 30, N0_c = 50),
    historical_control = list(N_t = 40, N0_c = 50))
  for (family in c("binomial", "normal")) {
    for (method in c("fixed", "mc")) {
      for (sizes in scenarios) {
        args <- sizes
        for (nm in names(sizes)) {
          suffix <- substring(nm, 2)
          if (family == "binomial") args[[paste0("y", suffix)]] <- sizes[[nm]] / 5 else {
            args[[paste0("mu", suffix)]] <- if (grepl("_t", nm)) 3 else 2
            args[[paste0("sigma", suffix)]] <- 1.5
          }
        }
        fit <- do.call(if (family == "binomial") bdpbinomial else bdpnormal,
          c(args, list(method = method, number_mcmc = 100, compare = FALSE)))
        component <- if (family == "binomial") "posterior" else "posterior_mu"
        draws <- posterior::as_draws_df(fit)
        expect_equal(draws$treatment, fit$posterior_treatment[[component]])
        if (fit$args1$arm2) {
          expect_equal(draws$control, fit$posterior_control[[component]])
          expect_equal(draws$difference, draws$treatment - draws$control)
        }
        expect_identical(posterior::variables(draws), tidy(fit)$term)
        expect_equal(nrow(draws), 100)
        expect_null(fit$final)
      }
    }
  }
})

test_that("draw formats retain regression coefficients, sigma, order, and actual draw counts", {
  skip_if_not_installed("posterior")
  for (fit in extractor_fits[c("lm_fixed", "lm_mc", "logit")]) {
    saved <- fit
    expected <- fit$posterior
    names(expected)[names(expected) == "intercept"] <- "(Intercept)"
    seed <- .Random.seed
    draws <- posterior::as_draws_df(fit)
    matrix_draws <- posterior::as_draws_matrix(fit)
    expect_s3_class(draws, "draws_df")
    expect_s3_class(matrix_draws, "draws_matrix")
    expect_equal(as.data.frame(draws)[names(expected)], expected, ignore_attr = TRUE)
    expect_identical(posterior::variables(matrix_draws), names(expected))
    expect_equal(dim(matrix_draws), dim(expected))
    expect_equal(as.numeric(matrix_draws), as.numeric(as.matrix(expected)))
    expect_equal(posterior::as_draws_df(matrix_draws), draws)
    expect_identical(draws$.chain, rep(1L, nrow(expected)))
    expect_identical(draws$.iteration, seq_len(nrow(expected)))
    expect_identical(draws$.draw, seq_len(nrow(expected)))
    expect_equal(posterior::ndraws(draws), glance(fit)$n.draws)
    expect_identical(.Random.seed, seed)
    expect_identical(fit, saved)
    expect_error(posterior::as_draws_df(fit, chains = 4), "Unsupported argument")
    expect_error(posterior::as_draws_matrix(fit, chains = 4), "Unsupported argument")
  }
})

test_that("survival draws preserve interval log hazards and reported summaries", {
  skip_if_not_installed("posterior")
  for (fit in extractor_fits[c("survival_one", "survival_two")]) {
    saved <- fit
    seed <- .Random.seed
    draws <- posterior::as_draws_df(fit)
    expect_equal(unname(as.matrix(as.data.frame(draws)[paste0("log_hazard_treatment[", 1:3, "]")])),
      unname(fit$posterior_treatment$posterior_log_hazard))
    if (fit$args1$arm2) {
      expect_equal(unname(as.matrix(as.data.frame(draws)[paste0("log_hazard_control[", 1:3, "]")])),
        unname(fit$posterior_control$posterior_log_hazard))
      expect_equal(mean(draws$treatment), tidy(fit)$estimate)
    } else {
      exposure <- c(2, 1, 0)
      expected <- exp(-drop(fit$posterior_treatment$posterior_hazard %*% exposure))
      expect_equal(draws$survival, expected)
      expect_equal(median(draws$survival), tidy(fit)$estimate)
      fit$args1$surv_time <- 0
      expect_equal(posterior::as_draws_df(fit)$survival, rep(1, nrow(draws)))
      fit <- saved
    }
    expect_identical(.Random.seed, seed)
    expect_identical(fit, saved)
    fit$posterior_treatment$posterior_log_hazard <- NULL
    if (fit$args1$arm2) fit$posterior_control$posterior_log_hazard <- NULL
    expect_equal(posterior::as_draws_df(fit), draws)
  }
})

test_that("modelsummary extracts posterior intervals and draws tables directly", {
  skip_if_not_installed("modelsummary")
  skip_if_not_installed("broom")
  withr::local_envvar(R_USER_CONFIG_DIR = tempdir())
  withr::local_seed(3052)
  fits <- c(list(
    binomial = bdpbinomial(y_t = 8, N_t = 40, y0_t = 12, N0_t = 60,
      method = "fixed", number_mcmc = 100),
    normal = bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
      mu0_t = 2.8, sigma0_t = 2, N0_t = 60,
      method = "fixed", number_mcmc = 100)), extractor_fits)
  for (fit in fits) {
    expected <- tidy(fit, conf.int = TRUE, conf.level = 0.8)
    actual <- modelsummary::get_estimates(fit, conf_level = 0.8)
    for (nm in c("term", "estimate", "std.error", "conf.low", "conf.high")) {
      expect_equal(actual[[nm]], expected[[nm]])
    }
    expect_true(all(is.na(actual$p.value)))
    expect_equal(tidy(fit, vcov = NULL, coef_rename = FALSE), tidy(fit))
    expect_equal(glance(fit, gof_map = "nobs"), glance(fit))
    expect_error(tidy(fit, vcov = diag(2)), "Covariance overrides")
    expect_error(tidy(fit, coef_rename = TRUE), "coef_rename must be FALSE")
    table <- modelsummary::modelsummary(fit, statistic = "conf.int", output = "data.frame")
    expect_s3_class(table, "data.frame")
    expect_true(any(table$statistic == "conf.int"))
    expect_true(any(table$term == "Num.Obs."))
  }
  table <- modelsummary::modelsummary(fits$lm_fixed, coef_map = c(treatment = "Treatment effect"),
    statistic = "conf.int", output = "data.frame")
  expect_true("Treatment effect" %in% table$term)
  for (fit in fits[c("logit", "survival_two")]) {
    actual <- modelsummary::get_estimates(fit, exponentiate = TRUE)
    expected <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    expect_equal(actual$estimate, expected$estimate)
    expect_equal(actual$conf.low, expected$conf.low)
    expect_equal(actual$conf.high, expected$conf.high)
    table <- modelsummary::modelsummary(fit, exponentiate = TRUE,
      statistic = "conf.int", output = "data.frame")
    expect_true(any(table$term == "Num.Obs."))
  }
})
