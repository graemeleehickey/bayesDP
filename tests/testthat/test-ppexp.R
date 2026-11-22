## Tests for the ppexp() piecewise-exponential CDF wrapper.

test_that("ppexp handles a single vector of hazard rates", {
  q <- 12
  x <- c(0.25, 0.3, 0.35, 0.4)
  cuts <- c(0, 6, 12, 18)

  pp <- ppexp(q, x, cuts)
  expect_length(pp, 1)
  expect_true(pp >= 0 && pp <= 1)
})

test_that("ppexp handles a matrix of hazard-rate draws", {
  set.seed(140)
  x <- matrix(rgamma(4 * 10, 0.1, 0.1), nrow = 10)
  cuts <- c(0, 6, 12, 18)

  pp <- ppexp(12, x, cuts)
  expect_length(pp, 10)
  expect_true(all(pp >= 0 & pp <= 1))
})

test_that("ppexp is monotonically non-decreasing in q", {
  x <- c(0.2, 0.3, 0.25, 0.4)
  cuts <- c(0, 5, 10, 15)
  qs <- c(1, 5, 10, 20)
  vals <- vapply(qs, ppexp, numeric(1), x = x, cuts = cuts)
  expect_true(all(diff(vals) >= 0))
})

test_that("ppexp rejects non-numeric, non-matrix input", {
  expect_error(
    ppexp(12, "a", c(0, 6, 12, 18)),
    "numeric vector or matrix"
  )
  expect_error(
    ppexp(12, list(1, 2), c(0, 6, 12, 18)),
    "numeric vector or matrix"
  )
})

test_that("ppexp_times matches per-time ppexp evaluation", {
  set.seed(8841)
  x <- matrix(rgamma(4 * 30, 0.1, 0.1), nrow = 30)
  cuts <- c(0, 6, 12, 18)
  times <- c(2, 5, 9, 14, 20)

  vectorized <- bayesDP:::ppexp_times(times, x, cuts)
  per_time <- lapply(times, ppexp, x, cuts = cuts)

  expect_length(vectorized, length(times))
  for (i in seq_along(times)) {
    # ppexpM() returns an n-by-1 column; ppexp_times() returns a plain vector.
    # The values must match (attributes such as dim are allowed to differ).
    expect_equal(vectorized[[i]], as.numeric(per_time[[i]]))
  }
})

test_that("piecewise-exponential probabilities agree with integrated hazards", {
  cuts <- c(0, 2, 5)
  hazards <- rbind(c(0.2, 0.4, 0.1), c(0.5, 0, 0.3))
  times <- c(0, 1, 2, 3, 5, 8)
  duration <- cbind(pmin(times, 2), pmax(0, pmin(times, 5) - 2),
    pmax(0, times - 5))
  expected <- -expm1(-duration %*% t(hazards))
  observed <- bayesDP:::ppexp_times(times, hazards, cuts)
  expect_equal(do.call(rbind, observed), unname(expected), tolerance = 1e-12)
  for (i in seq_along(times)) {
    expect_equal(as.numeric(ppexp(times[i], hazards, cuts)), expected[i, ])
    expect_equal(ppexp(times[i], hazards[1, ], cuts), expected[i, 1])
  }
})

test_that("a single interval gives the exponential CDF, including zero hazard", {
  for (rate in c(0, 0.2, 2)) {
    for (time in c(0, 1, 5)) {
      expect_equal(ppexp(time, rate, 0), -expm1(-rate * time), tolerance = 1e-12)
    }
  }
})
