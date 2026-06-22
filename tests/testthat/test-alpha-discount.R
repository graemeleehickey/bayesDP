## Tests for the educational helper alpha_discount().

test_that("alpha_discount matches the Weibull CDF discount functions", {
  # weibull: pweibull(p_hat, shape, scale) * alpha_max
  expect_equal(
    alpha_discount(0.1, discount_function = "weibull"),
    pweibull(0.1, shape = 3, scale = 0.135)
  )

  # scaledweibull: divides by pweibull(1, shape, scale)
  expect_equal(
    alpha_discount(0.1, discount_function = "scaledweibull"),
    pweibull(0.1, shape = 3, scale = 0.135) / pweibull(1, shape = 3, scale = 0.135)
  )

  # identity: returns p_hat directly
  expect_equal(alpha_discount(0.4, discount_function = "identity"), 0.4)
})

test_that("alpha_discount respects alpha_max and weibull parameters", {
  expect_equal(
    alpha_discount(0.2, discount_function = "weibull", alpha_max = 0.5),
    pweibull(0.2, shape = 3, scale = 0.135) * 0.5
  )

  expect_equal(
    alpha_discount(0.3,
      discount_function = "weibull",
      weibull_shape = 2, weibull_scale = 0.2
    ),
    pweibull(0.3, shape = 2, scale = 0.2)
  )

  # Output is a valid weight in [0, 1] for a default alpha_max.
  expect_true(alpha_discount(0.5) >= 0 && alpha_discount(0.5) <= 1)
})

test_that("alpha_discount rejects an unknown discount function", {
  expect_error(
    alpha_discount(0.5, discount_function = "not_a_function"),
    "discount_function input incorrectly"
  )
})
