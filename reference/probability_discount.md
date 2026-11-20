# Bayesian Discount Prior: Comparison Between Current and Historical Data

An educational calculation of the current-historical comparison for
normal summary data. This helper is not used internally. The fixed
method returns twice the smaller directional posterior probability; the
Monte Carlo method returns a normal-tail comparison for each draw.

## Usage

``` r
probability_discount(
  mu = NULL,
  sigma = NULL,
  N = NULL,
  mu0 = NULL,
  sigma0 = NULL,
  N0 = NULL,
  number_mcmc = 10000,
  method = "fixed"
)
```

## Arguments

- mu:

  scalar. Mean of the current data.

- sigma:

  scalar. Standard deviation of the current data.

- N:

  scalar. Number of observations of the current data.

- mu0:

  scalar. Mean of the historical data.

- sigma0:

  scalar. Standard deviation of the historical data.

- N0:

  scalar. Number of observations of the historical data.

- number_mcmc:

  scalar. Number of Monte Carlo simulations. Default is 10000.

- method:

  character. Default `"fixed"` estimates a directional probability q and
  returns `2 * min(q, 1 - q)`. With `"mc"`, each pair of mean draws is
  standardized by its conditional standard error and converted to a
  two-sided normal-tail comparison. See the normal vignette for details.

## Value

A numeric comparison in \[0, 1\]: one value under `"fixed"`, or a vector
of length `number_mcmc` under `"mc"`.

## Details

This is a borrowing comparison, not a posterior probability of equality
or equivalence. Weight variation under `"mc"` is part of that procedure;
more draws improve the precision of its distributional summaries. The
baseline prior for each source is proportional to \\1 / \sigma^2\\. The
fixed method uses direct Monte Carlo simulation.

## References

Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
Modeling and Simulation working group.(2017) Incorporation of stochastic
engineering models as prior information in Bayesian medical device
trials. *Journal of Biopharmaceutical Statistics*, 1-15.

## Examples

``` r
probability_discount(
  mu = 0, sigma = 1, N = 100,
  mu0 = 0.1, sigma0 = 1, N0 = 100
)
#> [1] 0.4832
```
