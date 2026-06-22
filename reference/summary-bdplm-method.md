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
function weight estimates are displayed as well. If `method`="mc", the
median estimate of alpha is displayed.
