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
the current data. Times are observed event or censoring times. Risk sets
and event counts describe current participants, not the split hazard
intervals. The survival estimates are posterior medians, with posterior
SDs and 95 equal-tailed credible intervals; they are not Kaplan–Meier
estimates. Saved fits without participant-level records must be refitted
before calling [`summary()`](https://rdrr.io/r/base/summary.html).

When a control arm is present, the summary describes the approximate
common log hazard ratio from the piecewise hazard draws. The comparison
between current and historical data and the resulting historical data
weight (alpha) are displayed, when historical data is present for the
respective arm. The displayed `coef` is the posterior mean log hazard
ratio, treatment relative to control; `exp(coef)` exponentiates this
mean, and `se(coef)` is its posterior SD. The 95 percent equal-tailed
credible limits are on the log hazard ratio scale. Sample sizes and
event counts describe the current arms, or historical controls when
current controls are absent.
