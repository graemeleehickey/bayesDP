test_that("survival response expressions give the same analysis as evaluated columns", {
  d <- data.frame(days = c(2, 4, 4, 8, 10, 12, 14, 16),
    event_code = c(2, 0, 2, 2, 0, 2, 0, 2), treatment = rep(0:1, 4))
  d$time <- d$days / 2
  d$status <- d$event_code == 2
  formulas <- list(
    Surv(days / 2, event_code == 2) ~ treatment,
    survival::Surv(days / 2, event_code == 2) ~ treatment,
    Surv(time = days / 2, event = event_code == 2) ~ treatment
  )
  for (method in c("fixed", "mc")) {
    fit_model <- function(formula) {
      set.seed(2815)
      bdpsurvival(formula, data = d, data0 = d, breaks = c(2, 5),
        method = method, number_mcmc = 500)
    }
    reference <- fit_model(Surv(time, status) ~ treatment)
    for (formula in formulas) {
      fit <- fit_model(formula)
      expect_identical(fit$posterior_treatment, reference$posterior_treatment)
      expect_identical(fit$posterior_control, reference$posterior_control)
      expect_identical(fit$args1$data, reference$args1$data)
      expect_identical(fit$args1$participants, reference$args1$participants)
      expect_equal(sum(fit$args1$data$exposure), 2 * sum(d$time))
      expect_equal(sum(fit$args1$data$status), 2 * sum(d$status))
      expect_equal(nrow(fit$args1$participants), 2 * nrow(d))
      expect_output(summary(fit), "Current treatment: n = 4")
    }
  }
  set.seed(2815)
  one <- bdpsurvival(survival::Surv(days / 2, event_code == 2) ~ 1,
    data = d, data0 = d, number_mcmc = 200)
  expect_equal(one$args1$participants$time, rep(d$time, 2))
})

test_that("invalid survival treatment values cannot silently exclude participants", {
  withr::local_seed(2816)
  rng <- .Random.seed
  d <- data.frame(time = 1:8, status = rep(0:1, 4), treatment = rep(0:1, 4))
  invalid <- list(c(0, 1, 2, 0, 1, 0, 1, 0), rep(1:2, 4),
    c(NA, d$treatment[-1]), c(Inf, d$treatment[-1]),
    as.character(d$treatment), factor(d$treatment), rep(.5, 8))
  for (source in c("Current", "Historical")) {
    for (treatment in invalid) {
      bad <- d; bad$treatment <- treatment
      expect_error(bdpsurvival(Surv(time, status) ~ treatment,
        data = if (source == "Current") bad else d,
        data0 = if (source == "Historical") bad else d), paste0(source, " treatment must contain only"))
      expect_identical(.Random.seed, rng)
    }
  }
  fit <- bdpsurvival(Surv(time, status) ~ treatment,
    data = transform(d, treatment = as.logical(treatment)), data0 = d, number_mcmc = 200)
  expect_equal(fit$args1$participants$treatment, rep(d$treatment, 2))
  expect_equal(sum(vapply(fit$args1[c("S_t", "S_c", "S0_t", "S0_c")],
    function(x) sum(x$status), numeric(1))), 2 * sum(d$status))
})

test_that("unsupported and incomplete survival responses stop without dropping rows", {
  d <- data.frame(time = 1:6, status = rep(0:1, 3), treatment = rep(0:1, 3), start = 0)
  for (formula in list(Surv(start, time, status) ~ treatment,
                       Surv(time, time + 1, type = "interval2") ~ treatment,
                       Surv(time, factor(status)) ~ treatment)) {
    expect_error(bdpsurvival(formula, data = d), "right-censored")
  }
  for (source in c("Current", "Historical")) {
    for (variable in c("time", "status")) {
      bad <- d; bad[[variable]][2] <- NA
      expect_error(bdpsurvival(Surv(time, status) ~ treatment,
        data = if (source == "Current") bad else d,
        data0 = if (source == "Historical") bad else d), paste0(source, " survival response"))
    }
  }
})
