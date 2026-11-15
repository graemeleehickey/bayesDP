# bdpbinomial Object Summary

`summary` method for class `bdpbinomial`.

## Usage

``` r
# S4 method for class 'bdpbinomial'
summary(object)
```

## Arguments

- object:

  object of class `bdpbinomial`. The result of a call to the
  [`bdpbinomial`](https://graemeleehickey.github.io/bayesDP/reference/bdpbinomial.md)
  function.

## Details

Displays input data, the current-historical comparison, and historical
weights when both sources are available. Under `method = "mc"`, the
comparison and weight summaries are medians across draws. The outcome
estimates are posterior medians of each arm's event probability.
Intervals are the 2.5th and 97.5th posterior percentiles. For a two-arm
analysis, the interval describes the treatment-minus-control difference,
calculated from paired draws. It is not a confidence interval or a test
of equality. Printed estimates retain four decimal places.
