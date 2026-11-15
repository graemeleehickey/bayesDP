# bdplm Object Summary

`summary` method for class `bdplm`.

## Usage

``` r
# S4 method for class 'bdplm'
summary(object)
```

## Arguments

- object:

  an object of class `bdplm`, a result of a call to
  [`bdplm`](https://graemeleehickey.github.io/bayesDP/reference/bdplm.md).

## Details

Displays a summary of the `bdplm` fit. Displayed output is similar to
[`summary.lm`](https://rdrr.io/r/stats/summary.lm.html). The function
call, residuals, and coefficient table are displayed. The discount
function weight estimates are displayed as well. If `method = "mc"`, the
median estimate of alpha is displayed. Residuals use current outcomes
minus predictions from posterior mean coefficients, with covariates and
intercept on the original scale. The method is also used by
[`bdplogit()`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md);
for those fits, residuals are observed binary outcomes minus fitted
probabilities obtained by applying the inverse logit to the posterior
mean linear predictor. Saved fits without the covariate centering means
must be refitted before calling
[`summary()`](https://rdrr.io/r/base/summary.html).
