# Count current participants

Return the number of current participants included in a fitted analysis.

## Usage

``` r
# S3 method for class 'bdpbinomial'
nobs(object, ...)

# S3 method for class 'bdpnormal'
nobs(object, ...)

# S3 method for class 'bdplm'
nobs(object, ...)

# S3 method for class 'bdpsurvival'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `bdpbinomial`, `bdpnormal`, `bdplm`, or `bdpsurvival` object.
  Fits from
  [`bdplogit()`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md)
  use class `bdplm`.

- ...:

  Additional arguments are not supported and produce an error.

## Value

A single numeric count, matching the `nobs` column of `glance(object)`.

## Details

Historical participants, including historical-only controls, are
excluded. These are unweighted participant counts, not effective sample
sizes. Survival counts use original participant records rather than
split follow-up intervals. For older survival or regression fits without
the necessary records, the result is `NA`. Use
[`glance()`](https://graemeleehickey.github.io/bayesDP/reference/glance.md)
to obtain historical counts separately.

## See also

[`glance`](https://graemeleehickey.github.io/bayesDP/reference/glance.md)

## Author

Graeme L. Hickey

## Examples

``` r
set.seed(42)
fit <- bdpbinomial(y_t = 8, N_t = 40, y0_t = 12, N0_t = 60,
  method = "fixed", number_mcmc = 100)
nobs(fit)
#> [1] 40
glance(fit)
#> # A tibble: 1 × 5
#>    nobs n.historical n.draws alpha.treatment alpha.control
#>   <dbl>        <dbl>   <dbl>           <dbl>         <dbl>
#> 1    40           60     100            0.82            NA
```
