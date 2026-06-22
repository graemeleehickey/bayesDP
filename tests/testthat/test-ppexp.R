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
