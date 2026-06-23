# Bayesian Discount Prior: Comparison Between Current and Historical Data

`probability_discount` can be used to estimate the posterior probability
of the comparison between historical and current data in the context of
a clinical trial with normal (mean) data. `probability_discount` is not
used internally but is given for educational purposes.

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

  character. Analysis method. Default value "`fixed`" estimates the
  posterior probability and holds it fixed. Alternative method "`mc`"
  estimates the posterior probability for each Monte Carlo iteration.
  See the the `bdpnormal` vignette  
  [`vignette("bdpnormal-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpnormal-vignette.md)
  for more details.

## Value

`probability_discount` returns an object of class
"probability_discount".

An object of class `probability_discount` contains the following:

- `p_hat`:

  scalar. The posterior probability of the comparison historical data
  weight. If `method="mc"`, a vector of posterior probabilities of
  length `number_mcmc` is returned.

## Details

This function is not used internally but is given for educational
purposes. Given the inputs, the output is the posterior probability of
the comparison between current and historical data in the context of a
clinical trial with normal (mean) data.

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
#> [1] 0.4858
```
