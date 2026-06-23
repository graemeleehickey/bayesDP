# Bayesian Discount Prior: Gaussian mean values

`bdpnormal` is used for estimating posterior samples from a Gaussian
outcome where an informative prior is used. The prior weight is
determined using a discount function. This code is modeled after the
methodologies developed in Haddad et al. (2017).

## Usage

``` r
bdpnormal(
  mu_t = NULL,
  sigma_t = NULL,
  N_t = NULL,
  mu0_t = NULL,
  sigma0_t = NULL,
  N0_t = NULL,
  mu_c = NULL,
  sigma_c = NULL,
  N_c = NULL,
  mu0_c = NULL,
  sigma0_c = NULL,
  N0_c = NULL,
  discount_function = "identity",
  alpha_max = 1,
  fix_alpha = FALSE,
  weibull_scale = 0.135,
  weibull_shape = 3,
  number_mcmc = 10000,
  method = "mc",
  compare = TRUE
)
```

## Arguments

- mu_t:

  scalar. Mean of the current treatment group.

- sigma_t:

  scalar. Standard deviation of the current treatment group.

- N_t:

  scalar. Number of observations of the current treatment group.

- mu0_t:

  scalar. Mean of the historical treatment group.

- sigma0_t:

  scalar. Standard deviation of the historical treatment group.

- N0_t:

  scalar. Number of observations of the historical treatment group.

- mu_c:

  scalar. Mean of the current control group.

- sigma_c:

  scalar. Standard deviation of the current control group.

- N_c:

  scalar. Number of observations of the current control group.

- mu0_c:

  scalar. Mean of the historical control group.

- sigma0_c:

  scalar. Standard deviation of the historical control group.

- N0_c:

  scalar. Number of observations of the historical control group.

- discount_function:

  character. Specify the discount function to use. Currently supports
  `weibull`, `scaledweibull`, and `identity`. The discount function
  `scaledweibull` scales the output of the Weibull CDF to have a max
  value of 1. The `identity` discount function uses the posterior
  probability directly as the discount weight. Default value is
  "`identity`".

- alpha_max:

  scalar. Maximum weight the discount function can apply. Default is 1.
  For a two-arm trial, users may specify a vector of two values where
  the first value is used to weight the historical treatment group and
  the second value is used to weight the historical control group.

- fix_alpha:

  logical. Fix alpha at alpha_max? Default value is FALSE.

- weibull_scale:

  scalar. Scale parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 0.135. For a two-arm trial, users may specify a vector of two
  values where the first value is used to estimate the weight of the
  historical treatment group and the second value is used to estimate
  the weight of the historical control group. Not used when
  `discount_function` = "identity".

- weibull_shape:

  scalar. Shape parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 3. For a two-arm trial, users may specify a vector of two
  values where the first value is used to estimate the weight of the
  historical treatment group and the second value is used to estimate
  the weight of the historical control group. Not used when
  `discount_function` = "identity".

- number_mcmc:

  scalar. Number of Monte Carlo simulations. Default is 10000.

- method:

  character. Analysis method with respect to estimation of the weight
  parameter alpha. Default method "`mc`" estimates alpha for each Monte
  Carlo iteration. Alternate value "`fixed`" estimates alpha once and
  holds it fixed throughout the analysis. See the the `bdpnormal`
  vignette  
  [`vignette("bdpnormal-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpnormal-vignette.md)
  for more details.

- compare:

  logical. Should a comparison object be included in the fit? For a
  one-arm analysis, the comparison object is simply the posterior chain
  of the treatment group parameter. For a two-arm analysis, the
  comparison object is the posterior chain of the treatment effect that
  compares treatment and control. If `compare=TRUE`, the comparison
  object is accessible in the `final` slot, else the `final` slot is
  `NULL`. Default is `TRUE`.

## Value

`bdpnormal` returns an object of class "bdpnormal". The functions
[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpnormal-method.md)
and
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpnormal-method.md)
are used to obtain and print a summary of the results, including user
inputs. The
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpnormal-method.md)
function displays visual outputs as well.

An object of class `bdpnormal` is a list containing at least the
following components:

- `posterior_treatment`:

  list. Entries contain values related to the treatment group:

  - `alpha_discount` numeric. Alpha value, the weighting parameter of
    the historical data.

  - `p_hat` numeric. The posterior probability of the stochastic
    comparison between the current and historical data.

  - `posterior_mu` vector. A vector of length `number_mcmc` containing
    the posterior mean of the treatment group. If historical treatment
    data is present, the posterior incorporates the weighted historical
    data.

  - `posterior_sigma2` vector. A vector of length `number_mcmc`
    containing the posterior variance of the treatment group. If
    historical treatment data is present, the posterior incorporates the
    weighted historical data.

  - `posterior_flat_mu` vector. A vector of length `number_mcmc`
    containing Monte Carlo samples of the mean of the current treatment
    group under a flat/non-informative prior, i.e., no incorporation of
    the historical data.

  - `posterior_flat_sigma2` vector. A vector of length `number_mcmc`
    containing Monte Carlo samples of the standard deviation of the
    current treatment group under a flat/non-informative prior, i.e., no
    incorporation of the historical data.

  - `prior_mu` vector. If historical treatment data is present, a vector
    of length `number_mcmc` containing Monte Carlo samples of the mean
    of the historical treatment group under a flat/non-informative
    prior.

  - `prior_sigma2` vector. If historical treatment data is present, a
    vector of length `number_mcmc` containing Monte Carlo samples of the
    standard deviation of the historical treatment group under a
    flat/non-informative prior.

- `posterior_control`:

  list. Similar entries as `posterior_treament`. Only present if a
  control group is specified.

- `final`:

  list. Contains the final comparison object, dependent on the analysis
  type:

  - One-arm analysis: vector. Posterior chain of the mean.

  - Two-arm analysis: vector. Posterior chain of the mean difference
    comparing treatment and control groups.

- `args1`:

  list. Entries contain user inputs. In addition, the following elements
  are output:

  - `arm2` binary indicator. Used internally to indicate one-arm or
    two-arm analysis.

  - `intent` character. Denotes current/historical status of treatment
    and control groups.

## Details

`bdpnormal` uses a two-stage approach for determining the strength of
historical data in estimation of a mean outcome. In the first stage, a
*discount function* is used that that defines the maximum strength of
the historical data and discounts based on disagreement with the current
data. Disagreement between current and historical data is determined by
stochastically comparing the respective posterior distributions under
noninformative priors. With Gaussian data, the comparison is the
proability (`p`) that the current mean is less than the historical mean.
The comparison metric `p` is then input into the discount function and
the final strength of the historical data is returned (alpha).

In the second stage, posterior estimation is performed where the
discount function parameter, `alpha`, is used incorporated in all
posterior estimation procedures.

To carry out a single arm (OPC) analysis, data for the current treatment
(`mu_t`, `sigma_t`, and `N_t`) and historical treatment (`mu0_t`,
`sigma0_t`, and `N0_t`) must be input. The results are then based on the
posterior distribution of the current data augmented by the historical
data.

To carry our a two-arm (RCT) analysis, data for the current treatment
and at least one of current or historical control data must be input.
The results are then based on the posterior distribution of the
difference between current treatment and control, augmented by available
historical data.

For more details, see the `bdpnormal` vignette:  
[`vignette("bdpnormal-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpnormal-vignette.md)

## References

Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
Modeling and Simulation working group.(2017) Incorporation of stochastic
engineering models as prior information in Bayesian medical device
trials. *Journal of Biopharmaceutical Statistics*, 1-15.

## See also

[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpnormal-method.md),
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpnormal-method.md),
and
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpnormal-method.md)
for details of each of the supported methods.

## Examples

``` r
# One-arm trial (OPC) example
fit <- bdpnormal(
  mu_t = 30, sigma_t = 10, N_t = 50,
  mu0_t = 32, sigma0_t = 10, N0_t = 50,
  method = "fixed"
)
summary(fit)
#> 
#>     One-armed bdp normal
#> 
#> data:
#>   Current treatment: mu_t = 30, sigma_t = 10, N_t = 50
#>   Historical treatment: mu0_t = 32, sigma0_t = 10, N0_t = 50
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.326
#> Discount function value (alpha) - treatment: 0.326
#> 95 percent CI: 
#>  27.9583  33.0088
#> posterior sample estimate:
#> mean of treatment group
#>  30.4735
if (FALSE) { # \dontrun{
plot(fit)
} # }

# Two-arm (RCT) example
fit2 <- bdpnormal(
  mu_t = 30, sigma_t = 10, N_t = 50,
  mu0_t = 32, sigma0_t = 10, N0_t = 50,
  mu_c = 25, sigma_c = 10, N_c = 50,
  mu0_c = 25, sigma0_c = 10, N0_c = 50,
  method = "fixed"
)
summary(fit2)
#> 
#>     Two-armed bdp normal
#> 
#> data:
#>   Current treatment: mu_t = 30, sigma_t = 10, N_t = 50
#>   Current control: mu_c = 25, sigma_c = 10, N_c = 50
#>   Historical treatment: mu0_t = 32, sigma0_t = 10, N0_t = 50
#>   Historical control: mu0_c = 25, sigma0_c = 10, N0_c = 50
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.33
#> Stochastic comparison (p_hat) - control (current vs. historical data): 0.993
#> Discount function value (alpha) - treatment: 0.33
#> Discount function value (alpha) - control: 0.993
#> alternative hypothesis: two.sided
#> 95 percent CI: 
#>  2.3734  8.6811
#> posterior sample estimates:
#> treatment group  control group
#>           30.50          25.00
if (FALSE) { # \dontrun{
plot(fit2)
} # }
```
