# bdpsurvival Object Plot

`plot` method for class `bdpsurvival`.

## Usage

``` r
# S4 method for class 'bdpsurvival'
plot(x, type = NULL, print = TRUE)
```

## Arguments

- x:

  object of class `bdpsurvival`. The result of a call to the
  [`bdpsurvival`](https://graemeleehickey.github.io/bayesDP/reference/bdpsurvival.md)
  function.

- type:

  character. Optional. Select plot type to print. Supports the
  following: "discount" gives the discount function and "survival" gives
  the survival curves. Leave NULL to display all plots in sequence.

- print:

  logical. Optional. Produce a plot (`print = TRUE`; default) or return
  a ggplot object (`print = FALSE`). When `print = FALSE`, it is
  possible to use `ggplot2` syntax to modify the plot appearance.

## Details

The `plot` method for `bdpsurvival` returns up to two plots: (1)
posterior survival curves and (2) the discount function. A call to
`plot` that omits the `type` input returns one plot at a time and
prompts the user to click the plot or press return to proceed to the
next plot. Otherwise, the user can specify a plot `type` to display the
requested plot.
