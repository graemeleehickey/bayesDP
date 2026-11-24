tidier_fits <- make_tidier_fits()

test_that("binomial and normal tidiers report paired contrasts in every input scenario", {
  withr::local_seed(2982)
  scenarios <- list(
    current_only = list(N_t = 40),
    one_arm = list(N_t = 40, N0_t = 60),
    two_arm = list(N_t = 40, N0_t = 60, N_c = 30, N0_c = 50),
    historical_control = list(N_t = 40, N0_c = 50)
  )
  for (model in c("binomial", "normal")) {
    for (method in c("fixed", "mc")) {
      for (sizes in scenarios) {
        args <- sizes
        for (nm in names(sizes)) {
          suffix <- substring(nm, 2)
          if (model == "binomial") {
            args[[paste0("y", suffix)]] <- sizes[[nm]] / 5
          } else {
            args[[paste0("mu", suffix)]] <- if (grepl("_t", nm)) 3 else 2
            args[[paste0("sigma", suffix)]] <- 1.5
          }
        }
        args <- c(args, list(method = method, number_mcmc = 300, compare = FALSE))
        fit <- do.call(if (model == "binomial") bdpbinomial else bdpnormal, args)
        component <- if (model == "binomial") "posterior" else "posterior_mu"
        t <- fit$posterior_treatment[[component]]
        c <- fit$posterior_control[[component]]
        draws <- if (is.null(c)) cbind(treatment = t) else cbind(treatment = t, control = c, difference = t - c)
        actual <- tidy(fit, conf.int = TRUE, conf.level = 0.8)
        expect_s3_class(actual, "tbl_df")
        expect_identical(actual$term, colnames(draws))
        expect_equal(actual$estimate, unname(apply(draws, 2, median)))
        expect_equal(actual$std.error, unname(apply(draws, 2, sd)))
        expect_equal(actual$conf.low, unname(apply(draws, 2, quantile, 0.1)))
        expect_equal(actual$conf.high, unname(apply(draws, 2, quantile, 0.9)))
        expect_named(tidy(fit), c("term", "estimate", "std.error"))
        g <- glance(fit)
        expect_named(g, c("nobs", "n.historical", "n.draws", "alpha.treatment", "alpha.control"))
        expect_equal(g$nobs, sum(sizes$N_t, sizes$N_c))
        expect_equal(g$n.historical, sum(sizes$N0_t, sizes$N0_c))
        expect_equal(g$n.draws, 300)
        for (arm in c("treatment", "control")) {
          alpha <- fit[[paste0("posterior_", arm)]]$alpha_discount
          expect_equal(g[[paste0("alpha.", arm)]], if (is.null(alpha)) NA_real_ else median(alpha))
        }
        expect_null(fit$final)
      }
    }
  }
})

test_that("linear and logistic regression summaries retain factors and posterior means", {
  for (fit in tidier_fits[c("lm_fixed", "lm_mc", "logit")]) {
    actual <- tidy(fit, conf.int = TRUE)
    expect_identical(actual$term, c("(Intercept)", "treatment", "x", "siteB"))
    draws <- fit$posterior[, c("intercept", "treatment", "x", "siteB")]
    expect_equal(actual$estimate, unname(colMeans(draws)))
    expect_equal(actual$std.error, unname(vapply(draws, sd, numeric(1))))
    expect_equal(actual$conf.low, unname(vapply(draws, quantile, numeric(1), 0.025)))
    expect_equal(actual$conf.high, unname(vapply(draws, quantile, numeric(1), 0.975)))
    g <- glance(fit)
    expect_equal(g$nobs, 48)
    expect_equal(g$n.historical, 48)
    expect_equal(g$n.draws, nrow(fit$posterior))
    expect_equal(g$alpha.treatment, median(fit$alpha_discount$treatment))
    expect_equal(g$alpha.control, median(fit$alpha_discount$control))
  }
  fit <- tidier_fits$lm_mc
  expect_gt(glance(fit)$n.draws, 301) # Actual retained size includes rounding up.
  sigma <- tidy(fit, effects = "sigma", conf.int = TRUE)
  expect_identical(sigma$term, "sigma")
  expect_equal(sigma$estimate, mean(fit$posterior$sigma))
  expect_equal(sigma$std.error, sd(fit$posterior$sigma))
  expect_equal(c(sigma$conf.low, sigma$conf.high), unname(quantile(fit$posterior$sigma, c(0.025, 0.975))))
  expect_error(tidy(tidier_fits$logit, effects = "sigma"), "no residual SD")
})

test_that("exponentiation has explicit point estimate, SD, and interval scales", {
  fit <- tidier_fits$logit
  log <- tidy(fit, conf.int = TRUE)
  ratio <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  expect_equal(ratio$estimate, exp(log$estimate))
  expect_equal(ratio$conf.low, exp(log$conf.low))
  expect_equal(ratio$conf.high, exp(log$conf.high))
  expect_equal(ratio$std.error, unname(vapply(fit$posterior, function(z) sd(exp(z)), numeric(1))))
  expect_false(isTRUE(all.equal(ratio$estimate, unname(colMeans(exp(fit$posterior))))))
  expect_error(tidy(tidier_fits$lm_fixed, exponentiate = TRUE), "only supported for logistic")
})

test_that("one-arm survival estimates use the requested horizon and original participants", {
  fit <- tidier_fits$survival_one
  for (time in c(0, 1, 2, 3, 5, 10)) {
    exposure <- c(min(time, 2), max(0, min(time, 5) - 2), max(0, time - 5))
    draws <- exp(-drop(fit$posterior_treatment$posterior_hazard %*% exposure))
    actual <- tidy(fit, time = time, conf.int = TRUE)
    expect_identical(actual$term, "survival")
    expect_equal(actual$time, time)
    expect_equal(actual$estimate, median(draws))
    expect_equal(actual$std.error, sd(draws))
    expect_equal(c(actual$conf.low, actual$conf.high), unname(quantile(draws, c(0.025, 0.975))))
  }
  expect_equal(tidy(fit), tidy(fit, time = 3))
  g <- glance(fit)
  expect_equal(g$nobs, 6)
  expect_equal(g$n.historical, 6)
  expect_equal(g$n.event, 4)
  expect_equal(g$n.event.historical, 4)
  expect_equal(g$n.draws, 300)
  expect_true(is.na(g$alpha.control))
  expect_gt(nrow(fit$args1$data), g$nobs + g$n.historical)
  expect_error(tidy(fit, exponentiate = TRUE), "only supported for two-arm")
})

test_that("two-arm survival uses the existing approximate effect even without final draws", {
  fit <- tidier_fits$survival_two
  differences <- fit$posterior_treatment$posterior_log_hazard - fit$posterior_control$posterior_log_hazard
  weights <- 1 / apply(differences, 2, var)
  draws <- drop(differences %*% (weights / sum(weights)))
  actual <- tidy(fit, conf.int = TRUE)
  expect_identical(actual$term, "treatment")
  expect_true(is.na(actual$time))
  expect_equal(actual$estimate, mean(draws))
  expect_equal(actual$std.error, sd(draws))
  expect_equal(c(actual$conf.low, actual$conf.high), unname(quantile(draws, c(0.025, 0.975))))
  ratio <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  expect_equal(ratio$estimate, exp(mean(draws)))
  expect_equal(ratio$std.error, sd(exp(draws)))
  expect_equal(c(ratio$conf.low, ratio$conf.high), exp(c(actual$conf.low, actual$conf.high)))
  expect_null(fit$final)
  expect_error(tidy(fit, time = 3), "only supported for one-arm")
  g <- glance(fit)
  expect_equal(g$nobs, 12)
  expect_equal(g$n.historical, 12)
  expect_equal(g$n.event, 8)
  expect_equal(g$n.event.historical, 8)
  expect_equal(g$alpha.treatment, median(fit$posterior_treatment$alpha_discount))
  expect_equal(g$alpha.control, median(fit$posterior_control$alpha_discount))
})

test_that("survival glance handles current-only and historical-only control scenarios", {
  withr::local_seed(2983)
  current <- data.frame(time = c(1, 2, 2, 4, 6, 8), status = c(1, 0, 1, 0, 1, 1), treatment = 1)
  historical <- data.frame(time = c(1, 3, 5), status = c(1, 0, 1), treatment = 0)
  expect_warning(one <- bdpsurvival(Surv(time, status) ~ 1, data = current,
    breaks = c(2, 5), method = "fixed", number_mcmc = 300), "Historical treatment data missing")
  expect_warning(two <- bdpsurvival(Surv(time, status) ~ treatment,
    data = current, data0 = historical, breaks = c(2, 5),
    method = "fixed", number_mcmc = 300), "Current control data missing")
  g1 <- glance(one)
  g2 <- glance(two)
  expect_identical(names(g1), names(g2))
  expect_equal(g1$n.historical, 0)
  expect_equal(g1$n.event.historical, 0)
  expect_equal(g2$nobs, 6)
  expect_equal(g2$n.historical, 3)
  expect_equal(g2$n.event, 4)
  expect_equal(g2$n.event.historical, 2)
  expect_true(all(is.na(g2[, c("alpha.treatment", "alpha.control")])))
  expect_equal(nrow(tidy(two)), 1L)
})

test_that("legacy fits retain available summaries without misleading counts", {
  for (fit in tidier_fits[c("survival_one", "survival_two")]) {
    original <- tidy(fit, conf.int = TRUE)
    fit$args1$participants <- NULL
    fit$posterior_treatment$posterior_log_hazard <- NULL
    if (!is.null(fit$posterior_control)) fit$posterior_control$posterior_log_hazard <- NULL
    expect_equal(tidy(fit, conf.int = TRUE), original)
    g <- glance(fit)
    expect_true(all(is.na(g[, c("nobs", "n.historical", "n.event", "n.event.historical")])))
    expect_equal(g$n.draws, 300)
  }
  fit <- tidier_fits$lm_fixed
  fit$args1$covariate_means <- NULL
  expect_equal(tidy(fit), tidy(tidier_fits$lm_fixed))
})

test_that("invalid and unused arguments fail explicitly", {
  fit <- tidier_fits$lm_fixed
  for (level in list(NA_real_, NaN, Inf, 0, 1, -0.1, "0.95", c(0.8, 0.9), NULL)) {
    expect_error(tidy(fit, conf.level = level), "conf.level must")
  }
  for (value in list(NA, 1, "TRUE", c(TRUE, FALSE), NULL)) {
    expect_error(tidy(fit, conf.int = value), "conf.int must")
    expect_error(tidy(fit, exponentiate = value), "exponentiate must")
  }
  for (time in list(NA_real_, Inf, -1, c(1, 2), "3", NULL)) {
    expect_error(tidy(tidier_fits$survival_one, time = time), "time must")
  }
  expect_error(tidy(fit, conf.lvel = 0.8), "Unsupported argument.*conf.lvel")
  expect_error(glance(fit, TRUE), "Unsupported argument.*unnamed")
  expect_error(tidy(fit, effects = "other"), "arg.*should be one of")
})

test_that("prescribed weights, saved draw counts, and random-number state are preserved", {
  withr::local_seed(2984)
  fit <- bdpbinomial(y_t = 8, N_t = 40, y0_t = 12, N0_t = 60,
    y_c = 6, N_c = 30, y0_c = 10, N0_c = 50,
    method = "fixed", fix_alpha = TRUE, alpha_max = c(0, 0.7), number_mcmc = 300)
  fit$args1$number_mcmc <- 999
  g <- glance(fit)
  expect_equal(g$n.draws, 300)
  expect_equal(g$alpha.treatment, 0)
  expect_equal(g$alpha.control, 0.7)
  for (x in c(list(fit), tidier_fits)) {
    saved <- x
    seed <- .Random.seed
    expect_equal(nrow(glance(x)), 1L)
    tidy(x, conf.int = TRUE)
    expect_identical(.Random.seed, seed)
    expect_identical(x, saved)
  }
})

test_that("broom and generics dispatch to the package methods", {
  skip_if_not_installed("broom")
  withr::local_seed(2985)
  bin <- bdpbinomial(y_t = 8, N_t = 40, method = "fixed", number_mcmc = 100)
  norm <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40, method = "fixed", number_mcmc = 100)
  for (fit in c(list(bin, norm), tidier_fits)) {
    expect_equal(broom::tidy(fit, conf.int = TRUE), tidy(fit, conf.int = TRUE))
    expect_equal(generics::tidy(fit), tidy(fit))
    expect_equal(broom::glance(fit), glance(fit))
    expect_equal(generics::glance(fit), glance(fit))
  }
})
