test_that("all supplied binomial and normal groups reject impossible summaries", {
  withr::local_seed(2817)
  rng <- .Random.seed
  for (normal in c(FALSE, TRUE)) {
    fit <- if (normal) bdpnormal else bdpbinomial
    base <- if (normal) list(mu_t = 3, sigma_t = 2, N_t = 20) else list(y_t = 3, N_t = 20)
    for (suffix in c("_t", "0_t", "_c", "0_c")) {
      valid <- if (normal) list(mu = 3, sigma = 2, N = 20) else list(y = 3, N = 20)
      invalid <- list(N = list(0, -1, 1.5, Inf, NA_real_, c(10, 20), matrix(20), "20", numeric()))
      if (normal) {
        invalid$N <- c(invalid$N, list(1))
        invalid$mu <- list(Inf, NaN, c(1, 2), matrix(1), "3")
        invalid$sigma <- list(0, -2, Inf, NA_real_, c(1, 2), matrix(2))
      } else {
        invalid$y <- list(-1, -.5, .5, 21, Inf, NA_real_, c(1, 2), matrix(1), "3")
      }
      for (field in names(invalid)) {
        for (value in invalid[[field]]) {
          group <- valid
          group[[field]] <- value
          names(group) <- paste0(names(group), suffix)
          args <- if (suffix == "_t") group else c(base, group)
          expect_error(do.call(fit, args), paste0(field, suffix))
          expect_identical(.Random.seed, rng)
        }
      }
    }
  }
})

test_that("summary-data borrowing settings are validated before drawing samples", {
  withr::local_seed(2818)
  rng <- .Random.seed
  invalid <- list(
    method = list("unknown", NA_character_, character(), c("fixed", "mc"), 1),
    discount_function = list("unknown", NA_character_, character(), c("identity", "weibull")),
    fix_alpha = list(NA, 1, logical(), c(TRUE, FALSE), matrix(TRUE)),
    compare = list(NA, 1, logical(), c(TRUE, FALSE), matrix(TRUE)),
    number_mcmc = list(0, 1, 2.5, Inf, NA_real_, c(10, 20), matrix(10)),
    alpha_max = list(-.1, 1.1, NA_real_, Inf, numeric(), c(.1, .2, .3), matrix(.5)),
    weibull_shape = list(0, -1, Inf, NA_real_, numeric(), c(1, 2, 3)),
    weibull_scale = list(0, -1, Inf, NA_real_, numeric(), c(1, 2, 3))
  )
  for (normal in c(FALSE, TRUE)) {
    fit <- if (normal) bdpnormal else bdpbinomial
    base <- if (normal) list(mu_t = 3, sigma_t = 2, N_t = 20) else list(y_t = 3, N_t = 20)
    for (field in names(invalid)) {
      for (value in invalid[[field]]) {
        args <- base
        if (field %in% c("weibull_shape", "weibull_scale")) args$discount_function <- "weibull"
        args[[field]] <- value
        expect_error(do.call(fit, args), field)
        expect_identical(.Random.seed, rng)
      }
    }
  }
  for (name in c("a0", "b0")) {
    for (value in list(0, -1, Inf, NA_real_, c(1, 2), matrix(1))) {
      args <- list(y_t = 3, N_t = 20); args[[name]] <- value
      expect_error(do.call(bdpbinomial, args), name)
    }
  }
})

test_that("valid boundaries, unused Weibull settings, and absent groups remain supported", {
  withr::local_seed(2819)
  for (events in 0:1) {
    for (method in c("fixed", "mc")) {
      fit <- bdpbinomial(y_t = events, N_t = 1, method = method, number_mcmc = 200,
        weibull_shape = NA_real_, weibull_scale = NA_real_)
      expect_true(all(is.finite(fit$posterior_treatment$posterior)))
      expect_null(fit$posterior_control)
    }
  }
  normal <- bdpnormal(mu_t = 0, sigma_t = 1, N_t = 2, number_mcmc = 200)
  expect_true(all(is.finite(normal$posterior_treatment$posterior_mu)))
  for (fit in list(bdpbinomial, bdpnormal)) {
    args <- if (identical(fit, bdpbinomial)) list(y_t = 3, N_t = 20, y0_c = 4, N0_c = 20) else
      list(mu_t = 3, sigma_t = 2, N_t = 20, mu0_c = 2, sigma0_c = 2, N0_c = 20)
    args <- c(args, list(discount_function = "weibull", fix_alpha = TRUE,
      alpha_max = c(0, 1), weibull_shape = NA_real_, weibull_scale = NA_real_, number_mcmc = 200))
    result <- do.call(fit, args)
    expect_true(result$args1$arm2)
    expect_length(result$final$posterior, 200)
  }
})
