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

Displays a summary of the `bdpbinomial` fit including the input data,
the stochastic comparison between current and historical data, and the
resulting historical data weight (alpha). If historical data is missing
then no stochastic comparison nor weight are displayed.

In the case of a one-arm analysis, the displayed 95 percent CI contains
the lower and upper limits of the augmented event rate of the current
data. The displayed `sample estimate` is the event rate of the current
data augmented by the historical data.

When a control arm is present, a two-arm analysis is carried out. Now,
the displayed 95 percent CI contains the lower and upper limits of the
difference between the treatment and control arms with the historical
data augmented to current data, if present. The displayed
`sample estimates` are the event rates of the treatment and control
arms, each of which are augmented when historical data are present.
