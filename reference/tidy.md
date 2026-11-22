# Tidy posterior summaries

Summarise saved posterior draws in a tibble, without refitting or
drawing new samples. These methods also work with
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html).

## Usage

``` r
tidy(x, ...)

# S3 method for class 'bdpbinomial'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'bdpnormal'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'bdplm'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  effects = c("coefficients", "sigma"),
  ...
)

# S3 method for class 'bdpsurvival'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  time = x$args1$surv_time,
  ...
)
```

## Arguments

- x:

  A fitted `bdpbinomial`, `bdpnormal`, `bdplm`, or `bdpsurvival` object.
  Fits from
  [`bdplogit()`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md)
  use class `bdplm`.

- ...:

  Additional arguments are not supported and produce an error.

- conf.int:

  Logical; include equal-tailed credible limits? Defaults to `FALSE`.

- conf.level:

  Credible interval probability, strictly between zero and one. Defaults
  to `0.95`.

- exponentiate:

  Logical; exponentiate logistic coefficients or the two-arm survival
  effect? Defaults to `FALSE`. Not supported for other quantities.

- effects:

  For regression fits, `"coefficients"` returns the intercept, treatment
  effect, and covariate coefficients; `"sigma"` returns the residual SD
  of a linear model. Logistic fits have no residual SD parameter.

- time:

  For a one-arm survival fit, a single finite, nonnegative time at which
  to summarise survival. Defaults to the fit's `surv_time`. Cannot be
  supplied for a two-arm fit.

## Value

A tibble with one row per quantity and columns `term`, `estimate`, and
`std.error`. Despite its conventional name, `std.error` is a posterior
SD, not a frequentist standard error or Monte Carlo standard error. If
`conf.int = TRUE`, `conf.low` and `conf.high` contain equal-tailed
credible limits, computed using
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) with its
default interpolation. Survival results also contain `time`. Results
retain full numerical precision. No test statistics or p-values are
reported.

## Details

Binomial and normal fits return `treatment`, then `control` and
`difference` when a control arm is present. These describe event
probabilities or population means, and the treatment-minus-control
difference. Estimates are posterior medians. Differences are calculated
from paired draws, including when controls are historical only or the
fit used `compare = FALSE`.

Regression fits return posterior means on the original covariate scale,
with the intercept named `(Intercept)`. Logistic coefficients are on the
log-odds scale unless `exponentiate = TRUE`.

One-arm survival fits return a `survival` row containing the posterior
median survival probability and a `time` column. Two-arm fits return a
`treatment` row for the approximate common log hazard ratio, treatment
relative to control, with `time = NA`. Its estimate is the posterior
mean of the inverse-variance-weighted interval log hazard ratio draws
used by [`summary()`](https://rdrr.io/r/base/summary.html); it is not a
time-specific hazard ratio or an assumption of proportional hazards.

With `exponentiate = TRUE`, the estimate is the exponential of the
posterior mean log effect, not the posterior mean odds ratio or hazard
ratio. Credible limits are exponentiated log-scale quantiles. The
`std.error` column is the SD of the exponentiated draws, so it is on the
reported ratio scale, not the log scale or a delta-method approximation.
For an exponentiated logistic intercept, the estimate describes baseline
odds.

## See also

[`glance`](https://graemeleehickey.github.io/bayesDP/reference/glance.md),
[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdplm-method.md)

## Author

Graeme L. Hickey

## Examples

``` r
set.seed(42)
fit <- bdpbinomial(y_t = 8, N_t = 40, y_c = 5, N_c = 40,
  y0_t = 10, N0_t = 50, y0_c = 8, N0_c = 50,
  method = "fixed", number_mcmc = 1000)
tidy(fit, conf.int = TRUE)
#> # A tibble: 3 × 5
#>   term       estimate std.error conf.low conf.high
#>   <chr>         <dbl>     <dbl>    <dbl>     <dbl>
#> 1 treatment    0.204     0.0425   0.127      0.293
#> 2 control      0.147     0.0415   0.0759     0.234
#> 3 difference   0.0555    0.0591  -0.0681     0.179
glance(fit)
#> # A tibble: 1 × 5
#>    nobs n.historical n.draws alpha.treatment alpha.control
#>   <dbl>        <dbl>   <dbl>           <dbl>         <dbl>
#> 1    80          100    1000           0.998         0.676
```
