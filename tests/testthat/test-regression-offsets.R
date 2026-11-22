offset_test_data <- function(binary = TRUE) {
  d <- expand.grid(treatment = 0:1, x = -1:1, z = -1:1, replicate = 1:4)
  d$Y <- if (binary) {
    d$replicate %% 2
  } else {
    2 + .7 * d$treatment + .4 * d$x + .2 * d$z +
      c(-1.1, -.3, .2, 1.2)[d$replicate]
  }
  d[, c("Y", "treatment", "x", "z")]
}

test_that("regressions reject offset calls before evaluating offsets or sampling", {
  withr::local_seed(2711)
  rng <- .Random.seed
  d <- offset_test_data()
  formulas <- list(
    Y ~ treatment + x + offset(2 * treatment),
    Y ~ treatment + x + offset(0 * z),
    Y ~ treatment + x + offset(0),
    Y ~ treatment + x + offset(z) + offset(-z),
    Y ~ treatment + x + stats::offset(z),
    Y ~ treatment + x + stats:::offset(z),
    Y ~ treatment + x + I(stats::offset(z)),
    Y ~ treatment + x + offset(rnorm(nrow(d)))
  )
  for (model in c("bdplm", "bdplogit")) {
    for (formula in formulas) {
      expect_error(get(model)(formula, data = d, data0 = d, method = "fixed",
        number_mcmc_alpha = stop("Comparison sampling must not start"),
        number_mcmc_beta = stop("Coefficient sampling must not start")),
        paste0(model, "() does not support offset() terms. No model has been fitted."),
        fixed = TRUE)
      expect_identical(.Random.seed, rng)
    }
  }
})

test_that("offset rejection applies to every supported method and borrowing setting", {
  d <- offset_test_data()
  priors <- list(prior_treatment_effect = 0, prior_control_effect = 0,
    prior_treatment_sd = 1, prior_control_sd = 1)
  for (model in c("bdplm", "bdplogit")) {
    for (method in if (model == "bdplm") c("fixed", "mc") else "fixed") {
      settings <- list(list(), list(alpha_max = 0), priors)
      if (method == "fixed") {
        settings <- c(settings,
          lapply(c(0, .4, 1), function(alpha) list(fix_alpha = TRUE, alpha_max = alpha)),
          list(c(priors, list(fix_alpha = TRUE, alpha_max = 0))))
      }
      for (setting in settings) {
        expect_error(do.call(get(model), c(list(
          formula = Y ~ treatment + x + offset(2 * treatment),
          data = d, data0 = d, method = method), setting)),
          paste0(model, "() does not support offset() terms. No model has been fitted."),
          fixed = TRUE)
      }
    }
  }
})

test_that("a covariate named offset remains an ordinary estimated covariate", {
  for (model in c("bdplm", "bdplogit")) {
    d <- offset_test_data(binary = model == "bdplogit")
    renamed <- d
    names(renamed)[names(renamed) == "z"] <- "offset"
    for (method in if (model == "bdplm") c("fixed", "mc") else "fixed") {
      settings <- list(method = method, number_mcmc_alpha = 300,
        number_mcmc_beta = 300)
      if (model == "bdplm") settings <- c(settings,
        list(number_mcmc_sigmagrid = 300, number_mcmc_sigma = 50))
      set.seed(2712)
      reference <- do.call(get(model), c(list(formula = Y ~ treatment + x + z,
        data = d, data0 = d), settings))
      set.seed(2712)
      observed <- do.call(get(model), c(list(formula = Y ~ treatment + x + offset,
        data = renamed, data0 = renamed), settings))
      names(observed$posterior)[names(observed$posterior) == "offset"] <- "z"
      expect_identical(observed$posterior, reference$posterior)
      expect_identical(observed$alpha_discount, reference$alpha_discount)
    }
  }
})
