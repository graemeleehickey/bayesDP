# Convert saved posterior draws

Use
[`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
or
[`posterior::as_draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)
to extract saved draws in a standard format for posterior calculations
and plotting. The optional posterior package must be installed.

## Usage

``` r
# S3 method for class 'bdpbinomial'
as_draws_matrix(x, ...)

# S3 method for class 'bdpnormal'
as_draws_matrix(x, ...)

# S3 method for class 'bdplm'
as_draws_matrix(x, ...)

# S3 method for class 'bdpsurvival'
as_draws_matrix(x, ...)

# S3 method for class 'bdpbinomial'
as_draws_df(x, ...)

# S3 method for class 'bdpnormal'
as_draws_df(x, ...)

# S3 method for class 'bdplm'
as_draws_df(x, ...)

# S3 method for class 'bdpsurvival'
as_draws_df(x, ...)
```

## Arguments

- x:

  A fitted `bdpbinomial`, `bdpnormal`, `bdplm`, or `bdpsurvival` object.
  Fits from
  [`bdplogit()`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md)
  use class `bdplm`.

- ...:

  Additional arguments are not supported and produce an error.

## Value

A
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
or
[`posterior::draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html)
object, with one row per retained draw. The data-frame format includes
`.chain`, `.iteration`, and `.draw` metadata columns.

## Details

Conversion preserves the order and number of saved draws and does not
refit, resample, or alter R's random-number state. Variable definitions
are:

- Binomial and normal fits:

  `treatment`, then `control` and `difference` when a control arm is
  present. These are event probabilities or population means and their
  paired treatment-minus-control differences, including historical-only
  controls and fits with `compare = FALSE`. The conditional mean
  variances stored by normal fits are not included: they are not
  residual-variance draws.

- Regression fits:

  The intercept `(Intercept)`, treatment effect, and covariate
  coefficients on the original covariate scale. Logistic coefficients
  are on the log-odds scale. Linear fits additionally include the
  residual SD `sigma`.

- Survival fits:

  Interval log hazards named `log_hazard_treatment[1]`,
  `log_hazard_treatment[2]`, and so on, with corresponding
  `log_hazard_control[j]` columns when a control arm is present.
  Interval order follows `x$args1$breaks`. One-arm fits also include
  `survival`, the survival probability at `x$args1$surv_time`; two-arm
  fits include `treatment`, the approximate common log hazard ratio
  reported by
  [`tidy()`](https://graemeleehickey.github.io/bayesDP/reference/tidy.md).
  This summary does not impose proportional hazards. Older fits without
  saved log hazards use the logarithms of their stored hazards.

Discount weights, comparison probabilities, and separate prior or
current-only draws are not included.

The draws are represented as one sequence (`.chain = 1`), with
consecutive `.iteration` and `.draw` indices. These indices describe
retained row order, not the original iteration numbers before burn-in or
thinning. Logistic fits contain a single retained MCMC sequence; the
other fits contain Monte Carlo draws from the package's sampling
calculations. Conversion does not create independent chains or establish
convergence. In particular, splitting these rows into artificial chains
is not a substitute for independent fits when assessing convergence.

## See also

[`tidy`](https://graemeleehickey.github.io/bayesDP/reference/tidy.md),
[`regression_extractors`](https://graemeleehickey.github.io/bayesDP/reference/regression_extractors.md)

## Author

Graeme L. Hickey

## Examples

``` r
if (requireNamespace("posterior", quietly = TRUE)) {
  set.seed(42)
  fit <- bdpbinomial(y_t = 8, N_t = 40, y_c = 5, N_c = 40,
    y0_t = 10, N0_t = 50, y0_c = 8, N0_c = 50,
    method = "fixed", number_mcmc = 100)
  draws <- posterior::as_draws_df(fit)
  mean(draws$difference > 0)
  posterior::as_draws_matrix(fit)
}
#> # A draws_matrix: 100 iterations, 1 chains, and 3 variables
#>     variable
#> draw treatment control difference
#>   1       0.23   0.124      0.111
#>   2       0.20   0.135      0.060
#>   3       0.18   0.120      0.056
#>   4       0.19   0.142      0.051
#>   5       0.26   0.109      0.152
#>   6       0.15   0.068      0.082
#>   7       0.18   0.139      0.045
#>   8       0.19   0.211     -0.017
#>   9       0.19   0.150      0.038
#>   10      0.21   0.146      0.060
#> # ... with 90 more draws
```
