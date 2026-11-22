test_that("model matrices preserve formula order with a reordered model frame", {
  d <- data.frame(Y = 1:6, x = c(2, 4, 1, 3, 6, 5), z = 6:1)
  formula <- terms(Y ~ z + x)
  frame <- model.frame(formula, d)
  reordered <- frame[c("x", "Y", "z")]
  attr(reordered, "terms") <- attr(frame, "terms")
  observed <- model.matrixBayes(formula, reordered, drop.baseline = TRUE)
  expect_equal(observed, model.matrix(formula, frame))

  names(reordered)[1] <- "unknown"
  expect_error(model.matrixBayes(formula, reordered, drop.baseline = TRUE),
    "model frame and formula mismatch", fixed = TRUE)
})

test_that("model matrices can obtain numeric variables from the formula environment", {
  Y <- 1:4
  x <- c(-2, -1, 1, 2)
  formula <- terms(Y ~ x)
  expect_equal(model.matrixBayes(formula), model.matrix(formula))
})

test_that("character and ordered covariates use the requested contrasts", {
  withr::local_options(contrasts = c("contr.treatment", "contr.poly"))
  d <- data.frame(Y = 1:6, site = rep(c("A", "B", "C"), 2))
  formula <- terms(Y ~ site)
  expect_warning(observed <- model.matrixBayes(formula, d, drop.baseline = TRUE),
    "variable 'site' converted to a factor", fixed = TRUE)
  expect_identical(colnames(observed), c("(Intercept)", "siteB", "siteC"))
  expect_equal(as.numeric(observed), as.numeric(cbind(1,
    rep(c(0, 1, 0), 2), rep(c(0, 0, 1), 2))))

  d$site <- ordered(d$site)
  ordered_matrix <- model.matrixBayes(formula, d, drop.baseline = TRUE)
  expect_identical(colnames(ordered_matrix), c("(Intercept)", "site.L", "site.Q"))
  expect_equal(as.numeric(ordered_matrix[, -1]),
    as.numeric(contr.poly(3)[rep(1:3, 2), ]))
})

test_that("explicit contrast matrices and functions are respected", {
  d <- data.frame(Y = 1:6, site = factor(rep(c("A", "B", "C"), 2)))
  formula <- terms(Y ~ site)
  custom <- contr.sum(3)
  by_matrix <- model.matrixBayes(formula, d, drop.baseline = TRUE,
    contrasts.arg = list(site = custom))
  by_name <- model.matrixBayes(formula, d, drop.baseline = TRUE,
    contrasts.arg = list(site = "contr.sum"))
  expect_equal(as.numeric(by_matrix), as.numeric(by_name))
  expect_identical(colnames(by_matrix), colnames(by_name))
  expect_identical(attr(by_name, "contrasts")$site, "contr.sum")
  expect_equal(as.numeric(by_matrix[, -1]), as.numeric(custom[rep(1:3, 2), ]))
  expect_equal(unname(attr(by_matrix, "contrasts")$site), unname(custom))

  expect_error(model.matrixBayes(formula, d, drop.baseline = TRUE,
    contrasts.arg = list(custom)), "invalid 'contrasts.arg' argument", fixed = TRUE)
  expect_warning(ignored <- model.matrixBayes(formula, d, drop.baseline = TRUE,
    contrasts.arg = list(absent = custom)), "absent, its contrast will be ignored")
  expect_equal(ignored, model.matrixBayes(formula, d, drop.baseline = TRUE))

  # Contrasts already attached to a factor take precedence over defaults.
  contrasts(d$site) <- custom
  expect_equal(model.matrixBayes(formula, d, drop.baseline = TRUE), by_matrix)
})

test_that("intercept-only model matrices retain the number of observations", {
  d <- data.frame(row.names = paste0("patient", 1:4))
  observed <- model.matrixBayes(terms(~ 1), d, drop.baseline = TRUE)
  expect_equal(dim(observed), c(4L, 1L))
  expect_equal(as.numeric(observed), rep(1, 4))
  expect_identical(colnames(observed), "(Intercept)")
})
