# Glance at a Bayesian discount prior fit

Return sample sizes, retained draw counts, and historical discount
weights in a one-row tibble. These methods also work with
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html).

## Usage

``` r
glance(x, ...)

# S3 method for class 'bdpbinomial'
glance(x, ...)

# S3 method for class 'bdpnormal'
glance(x, ...)

# S3 method for class 'bdplm'
glance(x, ...)

# S3 method for class 'bdpsurvival'
glance(x, ...)
```

## Arguments

- x:

  A fitted `bdpbinomial`, `bdpnormal`, `bdplm`, or `bdpsurvival` object.
  Fits from
  [`bdplogit()`](https://graemeleehickey.github.io/bayesDP/reference/bdplogit.md)
  use class `bdplm`.

- ...:

  Additional arguments are not supported and produce an error.

## Value

A one-row tibble containing:

- nobs:

  Number of current participants included in the analysis, summed across
  arms.

- n.historical:

  Number of historical participants included in the analysis, summed
  across arms. Historical-only controls are counted here, not in `nobs`.
  An absent source contributes zero. These are unweighted counts, not
  effective sample sizes.

- n.draws:

  Number of retained outcome posterior draws, counted from the saved
  draws rather than the requested simulation size. This is not a
  simulation effective sample size.

- alpha.treatment, alpha.control:

  Median stored historical discount weight for each arm. A prescribed or
  estimated fixed weight is returned unchanged. `NA` indicates that a
  weight is not defined, including historical-only controls and absent
  arms or historical sources.

- n.event, n.event.historical:

  For survival fits only, numbers of observed events among current and
  historical participants, respectively.

## Details

Survival sample sizes and event counts use the saved participant
records, not the rows created by splitting follow-up into hazard
intervals. For older survival fits without participant records, these
counts are `NA`; refit to obtain them. The columns and their order are
fixed within each fitted class. No likelihood criteria, convergence
diagnostics, or posterior comparison probabilities are supplied.

## See also

[`tidy`](https://graemeleehickey.github.io/bayesDP/reference/tidy.md)

## Author

Graeme L. Hickey

## Examples

``` r
set.seed(42)
fit <- bdpnormal(mu_t = 3, sigma_t = 2, N_t = 40,
  mu0_t = 2.8, sigma0_t = 2, N0_t = 60,
  method = "fixed", number_mcmc = 1000)
glance(fit)
#> # A tibble: 1 × 5
#>    nobs n.historical n.draws alpha.treatment alpha.control
#>   <dbl>        <dbl>   <dbl>           <dbl>         <dbl>
#> 1    40           60    1000           0.578            NA
```
