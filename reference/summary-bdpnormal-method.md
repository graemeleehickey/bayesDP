# bdpnormal Object Summary

`summary` method for class `bdpnormal`.

## Usage

``` r
# S4 method for class 'bdpnormal'
summary(object)
```

## Arguments

- object:

  object of class `bdpnormal`. The result of a call to the
  [`bdpnormal`](https://graemeleehickey.github.io/bayesDP/reference/bdpnormal.md)
  function.

## Details

Displays a summary of the `bdpnormal` fit including the input data, the
stochastic comparison between current and historical data, and the
resulting historical data weight (alpha). If historical data is missing
then no stochastic comparison nor weight are displayed.

In the case of a one-arm analysis, the displayed 95 percent CI contains
the lower and upper limits of the augmented mean value of the current
data. The displayed `mean of treatment group` is the mean of the current
data augmented by the historical data.

When a control arm is present, a two-arm analysis is carried out. Now,
the displayed 95 percent CI contains the lower and upper limits of the
difference between the treatment and control arms with the historical
data augmented to current data, if present. The displayed posterior
sample estimates are the mean of the treatment and control arms, each of
which are augmented when historical data are present.
