# bdpsurvival Object Summary

`summary` method for class `bdpsurvival`.

## Usage

``` r
# S4 method for class 'bdpsurvival'
summary(object)
```

## Arguments

- object:

  an object of class `bdpsurvival`, a result of a call to
  [`bdpsurvival`](https://graemeleehickey.github.io/bayesDP/reference/bdpsurvival.md).

## Details

Displays a summary of the `bdpsurvival` fit. The output is different,
conditional on a one- or two-arm survival analysis.

In the case of a one-arm analysis, the stochastic comparison between
current and historical data and the resulting historical data weight
(alpha) are displayed, followed by a survival table containing augmented
posterior estimates of the survival probability at each time point for
the current data.

When a control arm is present, a two-arm analysis is carried out. A
two-arm survival analysis is similar to a cox proportional hazards
analysis, and the displayed summary reflects this. First, the stochastic
comparison between current and historical data and the resulting
historical data weight (alpha) are displayed, when historical data is
present for the respective arm. The displayed `coef` value is the
log-hazard between the augmented treatment and control arms
(log(treatment) - log(control)). The lower and upper 95 percent interval
limits correspond to the `coef` value.
