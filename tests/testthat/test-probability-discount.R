## Tests for the educational helper probability_discount().

test_that("probability_discount returns a scalar probability for method = fixed", {
  set.seed(4471)
  p <- probability_discount(
    mu = 0, sigma = 1, N = 100,
    mu0 = 0.1, sigma0 = 1, N0 = 100
  )

  expect_length(p, 1)
  expect_true(p >= 0 && p <= 1)
})

test_that("probability_discount returns a vector for method = mc", {
  set.seed(8120)
  p <- probability_discount(
    mu = 0, sigma = 1, N = 100,
    mu0 = 0.1, sigma0 = 1, N0 = 100,
    method = "mc", number_mcmc = 2000
  )

  expect_length(p, 2000)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("probability_discount is near 1 when current and historical data agree", {
  # Identical inputs imply maximal agreement, so the transformed two-sided
  # probability should be close to 1.
  set.seed(2090)
  p <- probability_discount(
    mu = 5, sigma = 2, N = 500,
    mu0 = 5, sigma0 = 2, N0 = 500,
    number_mcmc = 1e5
  )
  expect_gt(p, 0.9)
})

test_that("probability_discount rejects an invalid method", {
  expect_error(
    probability_discount(
      mu = 0, sigma = 1, N = 10,
      mu0 = 0, sigma0 = 1, N0 = 10,
      method = "invalid"
    ),
    "method entered incorrectly"
  )
})
