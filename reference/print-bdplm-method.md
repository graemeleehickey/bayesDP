# bdplm Object Print

`print` method for class `bdplm`.

## Usage

``` r
# S4 method for class 'bdplm'
print(x)
```

## Arguments

- x:

  object of class `bdplm`. The result of a call to the
  [`bdplm`](https://graemeleehickey.github.io/bayesDP/reference/bdplm.md)
  function.

## Details

Displays a print of the `bdplm` fit and the initial function call. The
fit includes the estimate of the intercept, treatment effect, and
covariate effects. The discount function weight estimates are displayed
as well. If `method`="mc", the median estimate of alpha is displayed.
