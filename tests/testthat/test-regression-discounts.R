test_that("linear regression comparisons agree with the Student-t posterior", {
  withr::local_seed(2607)
  d <- data.frame(historical = rep(0:1, each = 60), x = rnorm(120))
  d$Y <- 2 + 0.4 * d$historical + 0.8 * d$x + rnorm(120, sd = 2)
  reference <- summary(lm(Y ~ historical + x, data = d))
  statistic <- reference$coefficients["historical", "t value"]
  expected <- 2 * pt(-abs(statistic), df = reference$df[2])
  observed <- discount_lm(d, discount_function = "identity", alpha_max = 0.7,
    fix_alpha = FALSE, number_mcmc_alpha = 30000,
    weibull_shape = 2, weibull_scale = 0.5, method = "fixed")
  expect_equal(observed$p_hat, expected, tolerance = 0.015, scale = 1)
  expect_equal(observed$alpha_discount, 0.7 * observed$p_hat)
})

test_that("regression discount functions transform the same comparison draws", {
  withr::local_seed(2608)
  d <- data.frame(historical = rep(0:1, each = 60), x = rnorm(120))
  continuous <- transform(d, Y = 2 + 0.3 * historical + x + rnorm(120, sd = 2))
  binary <- transform(d, Y = rbinom(120, 1, plogis(-0.4 + 0.3 * historical + 0.2 * x)))
  models <- list(
    list(fit = discount_lm, data = continuous, methods = c("fixed", "mc")),
    list(fit = discount_logit, data = binary, methods = "fixed")
  )
  for (model in models) {
    for (method in model$methods) {
      calculate <- function(discount) {
        set.seed(2609)
        model$fit(model$data, discount_function = discount, alpha_max = 0.65,
          fix_alpha = FALSE, number_mcmc_alpha = 1000,
          weibull_shape = 2, weibull_scale = 0.5, method = method)
      }
      identity <- calculate("identity")
      expect_length(identity$p_hat, if (method == "mc") 1000 else 1)
      expect_true(all(identity$p_hat >= 0 & identity$p_hat <= 1))
      expect_equal(identity$alpha_discount, 0.65 * identity$p_hat)
      for (discount in c("weibull", "scaledweibull")) {
        observed <- calculate(discount)
        expected <- 1 - exp(-(identity$p_hat / 0.5)^2)
        if (discount == "scaledweibull") expected <- expected / (1 - exp(-4))
        expect_equal(observed$p_hat, identity$p_hat)
        expect_equal(observed$alpha_discount, 0.65 * expected)
      }
    }
  }
})

test_that("regression fits honour different prescribed weights for the two arms", {
  withr::local_seed(2610)
  make_data <- function(n, binary) {
    treatment <- rep(0:1, each = n / 2)
    x <- rnorm(n)
    eta <- -0.5 + 0.6 * treatment + 0.2 * x
    Y <- if (binary) rbinom(n, 1, plogis(eta)) else eta + rnorm(n)
    data.frame(Y = Y, treatment = treatment, x = x)
  }
  for (binary in c(FALSE, TRUE)) {
    current <- make_data(80, binary)
    historical <- make_data(120, binary)
    args <- list(formula = Y ~ treatment + x, data = current, data0 = historical,
      method = "fixed", fix_alpha = TRUE, alpha_max = c(0.6, 0.25),
      number_mcmc_alpha = 300, number_mcmc_beta = 300)
    if (!binary) args <- c(args,
      list(number_mcmc_sigmagrid = 300, number_mcmc_sigma = 50))
    fit <- do.call(if (binary) bdplogit else bdplm, args)
    expect_equal(fit$alpha_discount, data.frame(treatment = 0.6, control = 0.25))
    expect_equal(nrow(fit$posterior), 300)
    expect_true(all(is.finite(as.matrix(fit$posterior))))
    expect_named(fit$estimates$coefficients,
      c("intercept", "treatment", "x", if (!binary) "sigma"))
  }
})
