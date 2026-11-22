test_that("normal and binomial plots return and display each supported view", {
  withr::local_seed(2604)
  factories <- list(
    function(history) bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
      mu0_t = if (history) 3.2 else NULL,
      sigma0_t = if (history) 2 else NULL, N0_t = if (history) 50 else NULL,
      mu_c = 2, sigma_c = 2, N_c = 40, method = "fixed", number_mcmc = 300),
    function(history) bdpbinomial(y_t = 8, N_t = 40,
      y0_t = if (history) 10 else NULL, N0_t = if (history) 50 else NULL,
      y_c = 12, N_c = 40, method = "fixed", number_mcmc = 300)
  )
  for (factory in factories) {
    for (history in c(FALSE, TRUE)) {
      fit <- factory(history)
      views <- plot(fit, print = FALSE)
      expected <- c("posteriors", "density", if (history) "discount")
      expect_named(views, expected)
      expect_setequal(as.character(views$density$data$group), c("Treatment", "Control"))
      expect_equal(unique(as.character(views$density$data$information_sources)),
        "Augmented posterior")
      for (type in expected) {
        expect_equal(plot(fit, type = type, print = FALSE)$data, views[[type]]$data)
        expect_error(plot(fit, type = type), NA)
      }
      # The full sequence asks for input in interactive R sessions.
      if (!interactive()) {
        before <- par("ask")
        expect_error(plot(fit), NA)
        expect_identical(par("ask"), before)
      }
      if (!history) {
        expect_null(plot(fit, type = "discount", print = FALSE))
        expect_null(plot(fit, type = "discount"))
      }
      effect <- plot(fit, type = "effect", print = FALSE)
      expect_equal(plot(fit, type = "effect")$data, effect$data)
    }
  }
  one_arm <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40, number_mcmc = 300)
  expect_error(plot(one_arm, type = "effect"), "requires treatment and control data")
})

test_that("survival plots display the available sources and restore graphics settings", {
  withr::local_seed(2605)
  d <- data.frame(time = 1:6, status = c(1, 0, 1, 1, 0, 1))
  for (history in c(FALSE, TRUE)) {
    fit_model <- function() bdpsurvival(Surv(time, status) ~ 1,
      data = d, data0 = if (history) d else NULL,
      breaks = c(2, 4), number_mcmc = 300)
    if (history) {
      fit <- fit_model()
    } else {
      expect_warning(fit <- fit_model(), "Historical treatment data missing")
    }
    views <- plot(fit, print = FALSE)
    expect_named(views, c("survival", if (history) "discount"))
    expected_sources <- c("Current-only posterior", "Augmented posterior",
      if (history) "Historical-only posterior")
    expect_setequal(views$survival$data$source, expected_sources)
    expect_true(all(views$survival$data$y >= 0 & views$survival$data$y <= 1))
    for (type in names(views)) {
      expect_equal(plot(fit, type = type, print = FALSE)$data, views[[type]]$data)
      expect_error(plot(fit, type = type), NA)
    }
    if (!interactive()) {
      before <- par("ask")
      expect_error(plot(fit), NA)
      expect_identical(par("ask"), before)
    }
    if (!history) {
      expect_null(plot(fit, type = "discount", print = FALSE))
      expect_null(plot(fit, type = "discount"))
    }
  }
})

test_that("survival curves use the available control source without artificial borrowing", {
  withr::local_seed(2606)
  d <- data.frame(time = 1:6, status = c(1, 0, 1, 1, 0, 1), treatment = 1)
  control <- transform(d, time = time + 1, treatment = 0)
  for (historical_control in c(FALSE, TRUE)) {
    current <- if (historical_control) d else rbind(d, control)
    historical <- if (historical_control) rbind(d, control) else d
    fit_model <- function() bdpsurvival(Surv(time, status) ~ treatment,
      data = current, data0 = historical, breaks = c(2, 4), number_mcmc = 300)
    if (historical_control) {
      expect_warning(fit <- fit_model(), "Current control data missing")
    } else {
      fit <- fit_model()
    }
    curves <- plot(fit, type = "survival", print = FALSE)$data
    curves <- curves[curves$group == "Control", ]
    source <- if (historical_control) "Historical-only posterior" else "Current-only posterior"
    expect_setequal(curves$source, c(source, "Augmented posterior"))
    original <- curves[curves$source == source, c("x", "y")]
    augmented <- curves[curves$source == "Augmented posterior", c("x", "y")]
    rownames(original) <- rownames(augmented) <- NULL
    expect_equal(augmented, original)
    # Integrate the interval hazards independently to check the plotted values.
    lengths <- cbind(pmin(augmented$x, 2),
      pmax(0, pmin(augmented$x, 4) - 2), pmax(0, augmented$x - 4))
    expected <- apply(exp(-fit$posterior_control$posterior_hazard %*% t(lengths)),
      2, median)
    expect_equal(augmented$y, as.numeric(expected), tolerance = 1e-12)
  }
})
