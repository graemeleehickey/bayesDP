regression_fit_data <- function(binary) {
  d <- expand.grid(treatment = 0:1, x = -1:1, replicate = 1:6)
  d$Y <- if (binary) d$replicate %% 2 else
    2 + .7 * d$treatment + .3 * d$x + c(-1.2, -.6, -.1, .2, .7, 1)[d$replicate]
  d[, c("Y", "treatment", "x")]
}

test_that("treatment-only regressions return exactly the intended parameters", {
  withr::local_seed(2812)
  for (model in c("bdplm", "bdplogit")) {
    d <- regression_fit_data(model == "bdplogit")
    for (method in if (model == "bdplm") c("fixed", "mc") else "fixed") {
      for (prescribed in if (method == "fixed") c(FALSE, TRUE) else FALSE) {
        args <- list(formula = Y ~ treatment, data = d, data0 = d,
          method = method, fix_alpha = prescribed, alpha_max = c(.6, .3),
          number_mcmc_alpha = 500, number_mcmc_beta = 500)
        if (model == "bdplm") args <- c(args,
          list(number_mcmc_sigmagrid = 500, number_mcmc_sigma = 100))
        fit <- do.call(get(model), args)
        expected <- c("intercept", "treatment", if (model == "bdplm") "sigma")
        expect_named(fit$posterior, expected)
        expect_named(fit$estimates$coefficients, expected)
        expect_equal(nrow(fit$posterior), 500)
        expect_true(all(is.finite(as.matrix(fit$posterior))))
        expect_output(summary(fit), "Coefficients:")
        expect_output(print(fit), "Coefficients:")
        expect_equal(generics::tidy(fit)$term, c("(Intercept)", "treatment"))
      }
    }
  }
})

test_that("logistic fits use reproducible and distinct stage seeds from R", {
  withr::local_seed(2813)
  d <- regression_fit_data(TRUE)
  fit_model <- function(...) bdplogit(Y ~ treatment + x, data = d, data0 = d,
    number_mcmc_alpha = 500, number_mcmc_beta = 500, ...)
  before <- .Random.seed
  a <- fit_model()
  after <- .Random.seed
  expect_false(identical(before, after))
  expect_named(a$args1$mcmc_seeds, c("comparison_treatment", "comparison_control", "final"))
  expect_equal(length(unique(a$args1$mcmc_seeds)), 3)
  set.seed(2813)
  b <- fit_model(method = "fixed")
  expect_identical(a$posterior, b$posterior)
  expect_identical(a$alpha_discount, b$alpha_discount)
  expect_identical(a$args1$mcmc_seeds, b$args1$mcmc_seeds)
  expect_identical(after, .Random.seed)
  c <- fit_model()
  expect_false(identical(b$posterior, c$posterior))
  expect_false(identical(b$args1$mcmc_seeds, c$args1$mcmc_seeds))
  set.seed(2814)
  different <- fit_model()
  expect_false(identical(a$posterior, different$posterior))
  set.seed(2813)
  prescribed <- fit_model(fix_alpha = TRUE, alpha_max = .4)
  expect_identical(a$args1$mcmc_seeds, prescribed$args1$mcmc_seeds)
  expect_identical(after, .Random.seed)
  set.seed(2813)
  expect_identical(prescribed$posterior, fit_model(fix_alpha = TRUE, alpha_max = .4)$posterior)
})

test_that("logistic method validation rejects unsupported and malformed choices", {
  d <- regression_fit_data(TRUE)
  for (method in list("mc", "unknown", NA_character_, NULL, character(), c("fixed", "mc"), 1)) {
    expect_error(bdplogit(Y ~ treatment, data = d, data0 = d, method = method),
      "[Mm]ethod|method")
  }
})
