# bdpsurvival Object Print

`print` method for class `bdpsurvival`.

## Usage

``` r
# S4 method for class 'bdpsurvival'
print(x)
```

## Arguments

- x:

  object of class `bdpsurvival`. The result of a call to the
  [`bdpsurvival`](https://graemeleehickey.github.io/bayesDP/reference/bdpsurvival.md)
  function.

## Details

Displays a print of the `bdpsurvival` fit. The output is different,
conditional on a one- or two-arm survival analysis.

In the case of a one-arm analysis, a brief summary is displayed,
including the current data sample size, number of events, user input
survival time, the augmented median survival probability, and
corresponding lower and upper 95 percent interval limits.

When a control arm is present, the output is the same as a call to
[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpsurvival-method.md).
