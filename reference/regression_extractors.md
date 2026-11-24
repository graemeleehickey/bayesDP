# Extract regression estimates and credible intervals

Extract posterior summaries from linear and logistic regression fits
without refitting or drawing new samples.

## Usage

``` r
# S3 method for class 'bdplm'
coef(object, ...)

# S3 method for class 'bdplm'
vcov(object, complete = TRUE, ...)

# S3 method for class 'bdplm'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fit from
  [`bdplm`](https://graemeleehickey.github.io/bayesDP/reference/bdplm.md)
  or
  [`bdplogit`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md).
  Both use class `bdplm`.

- ...:

  Additional arguments are not supported and produce an error.

- complete:

  Logical, for compatibility with
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html). Both values
  return the same matrix because these fits do not contain aliased
  coefficients.

- parm:

  Coefficient names or positive integer positions. Defaults to all
  coefficients.

- level:

  Credible interval probability, strictly between zero and one.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector; [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns a square
matrix with matching coefficient names;
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per selected coefficient and two columns giving the lower
and upper posterior quantiles.

## Details

Coefficients are posterior means on the original covariate scale, with
the intercept named `(Intercept)`. Logistic coefficients remain on the
log-odds scale. The residual SD `sigma` is excluded; use
`tidy(object, effects = "sigma")` for that parameter in a linear model.

Despite the conventional function names,
[`confint()`](https://rdrr.io/r/stats/confint.html) returns equal-tailed
posterior credible intervals, using
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) with its
default interpolation. [`vcov()`](https://rdrr.io/r/stats/vcov.html)
returns the sample covariance of the saved coefficient draws, including
their posterior correlations. It is neither a frequentist sampling
covariance nor the Monte Carlo covariance of the estimated posterior
means.

## See also

[`tidy`](https://graemeleehickey.github.io/bayesDP/reference/tidy.md),
[`posterior_draws`](https://graemeleehickey.github.io/bayesDP/reference/posterior_draws.md),
[`nobs.bdpbinomial`](https://graemeleehickey.github.io/bayesDP/reference/nobs.bdpbinomial.md)

## Author

Graeme L. Hickey

## Examples

``` r
set.seed(42)
current <- data.frame(treatment = rep(0:1, each = 20))
current$y <- 2 + current$treatment + rnorm(40)
historical <- current
historical$y <- 2 + historical$treatment + rnorm(40)
fit <- bdplm(y ~ treatment, data = current, data0 = historical,
  method = "fixed", number_mcmc_alpha = 100,
  number_mcmc_sigmagrid = 100, number_mcmc_sigma = 50,
  number_mcmc_beta = 200)
coef(fit)
#> (Intercept)   treatment 
#>   2.0827151   0.7207231 
vcov(fit)
#>             (Intercept)   treatment
#> (Intercept)  0.03975653 -0.04218306
#> treatment   -0.04218306  0.10341119
confint(fit, parm = "treatment", level = 0.9)
#>                5%      95%
#> treatment 0.16551 1.231625
```
