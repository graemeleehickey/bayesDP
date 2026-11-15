# bdpnormal Object Plot

`plot` method for class `bdpnormal`.

## Usage

``` r
# S4 method for class 'bdpnormal'
plot(x, type = NULL, print = TRUE)
```

## Arguments

- x:

  object of class `bdpnormal`. The result of a call to the
  [`bdpnormal`](https://graemeleehickey.github.io/bayesDP/reference/bdpnormal.md)
  function.

- type:

  character. Optional. Select plot type to print. Supports the
  following: "discount" gives the discount function; "posteriors" gives
  the posterior distributions by information source; "density" gives the
  augmented arm distributions; and "effect" gives the
  treatment-minus-control difference for a two-arm analysis. Leave NULL
  to display the first three plots in sequence.

- print:

  logical. Optional. Produce a plot (`print = TRUE`; default) or return
  a ggplot object (`print = FALSE`). When `print = FALSE`, it is
  possible to use `ggplot2` syntax to modify the plot appearance.

## Details

The `plot` method for `bdpnormal` returns up to three plots: (1)
posterior density curves by information source; (2) augmented arm
densities; and (3) the discount function. A call to `plot` that omits
the `type` input returns one plot at a time and prompts the user to
click the plot or press return to proceed to the next plot. Otherwise,
the user can specify a plot `type` to display the requested plot.
