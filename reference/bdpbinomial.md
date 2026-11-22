# Bayesian Discount Prior: Binomial counts

`bdpbinomial` is used for estimating posterior samples from a binomial
outcome where an informative prior is used. The prior weight is
determined using a discount function. This code is modeled after the
methodologies developed in Haddad et al. (2017).

## Usage

``` r
bdpbinomial(
  y_t = NULL,
  N_t = NULL,
  y0_t = NULL,
  N0_t = NULL,
  y_c = NULL,
  N_c = NULL,
  y0_c = NULL,
  N0_c = NULL,
  discount_function = "identity",
  alpha_max = 1,
  fix_alpha = FALSE,
  a0 = 1,
  b0 = 1,
  number_mcmc = 10000,
  weibull_scale = 0.135,
  weibull_shape = 3,
  method = "mc",
  compare = TRUE
)
```

## Arguments

- y_t:

  scalar. Number of events for the current treatment group.

- N_t:

  scalar. Sample size of the current treatment group.

- y0_t:

  scalar. Number of events for the historical treatment group.

- N0_t:

  scalar. Sample size of the historical treatment group.

- y_c:

  scalar. Number of events for the current control group.

- N_c:

  scalar. Sample size of the current control group.

- y0_c:

  scalar. Number of events for the historical control group.

- N0_c:

  scalar. Sample size of the historical control group.

- discount_function:

  character. Specify the discount function to use. Currently supports
  `weibull`, `scaledweibull`, and `identity`. The discount function
  `scaledweibull` scales the output of the Weibull CDF to have a max
  value of 1. The `identity` discount function uses the transformed
  comparison directly as the discount weight. Default value is
  "`identity`".

- alpha_max:

  scalar in \[0, 1\]. Maximum weight the discount function can apply.
  Default is 1. For a two-arm trial, users may specify a vector of two
  values where the first value is used to weight the historical
  treatment group and the second value is used to weight the historical
  control group.

- fix_alpha:

  logical. Fix alpha at alpha_max? Default value is FALSE.

- a0:

  positive, finite scalar. Prior value for a beta shape parameter.
  Default is 1.

- b0:

  positive, finite scalar. Prior value for a beta shape parameter.
  Default is 1.

- number_mcmc:

  integer of at least 2. Number of Monte Carlo simulations. Default is
  10000.

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

- method:

  character. Analysis method with respect to estimation of the weight
  parameter alpha. Default method "`mc`" estimates alpha for each Monte
  Carlo iteration, giving a distribution of generated weights.
  Increasing the number of draws improves numerical precision but does
  not remove the spread of that distribution. These weights are not
  draws from a separately specified prior/posterior model for alpha.
  Alternate value "`fixed`" estimates alpha once and holds it fixed
  throughout the analysis. See the `bdpbinomial` vignette  
  [`vignette("bdpbinomial-vignette", package = "bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpbinomial-vignette.md)
  for more details.

- compare:

  logical. Should a comparison object be included in the fit? For a
  one-arm analysis, the comparison object is simply the posterior draws
  of the treatment group parameter. For a two-arm analysis, the
  comparison object is the posterior draws of the treatment effect that
  compares treatment and control. If `compare = TRUE`, the comparison
  object is accessible in the `final` slot, else the `final` slot is
  `NULL`. Default is `TRUE`.

## Value

`bdpbinomial` returns an object of class "bdpbinomial". The functions
[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpbinomial-method.md)
and
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpbinomial-method.md)
are used to obtain and print a summary of the results, including user
inputs. The
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpbinomial-method.md)
function displays visual outputs as well.

An object of class `bdpbinomial` is a list containing at least the
following components:

- `posterior_treatment`:

  list. Entries contain values related to the treatment group:

  - `alpha_discount` numeric. Alpha value, the weighting parameter of
    the historical data.

  - `p_hat` numeric. The current-historical comparison between the
    current and historical data.

  - `posterior` vector. A vector of length `number_mcmc` containing
    posterior Monte Carlo samples of the event rate of the treatment
    group. If historical treatment data is present, the posterior
    incorporates the weighted historical data.

  - `posterior_flat` vector. A vector of length `number_mcmc` containing
    Monte Carlo samples of the event rate of the current treatment group
    under a flat/non-informative prior, i.e., no incorporation of the
    historical data.

  - `prior` vector. If historical treatment data is present, a vector of
    length `number_mcmc` containing Monte Carlo samples of the event
    rate of the historical treatment group under a flat/non-informative
    prior.

- `posterior_control`:

  list. Similar entries as `posterior_treament`. Only present if a
  control group is specified.

- `final`:

  list. Contains the final comparison object, dependent on the analysis
  type:

  - One-arm analysis: vector. Posterior chain of binomial proportion.

  - Two-arm analysis: vector. Posterior chain of binomial proportion
    difference comparing treatment and control groups.

- `args1`:

  list. Entries contain user inputs. In addition, the following elements
  are output:

  - `arm2` binary indicator. Used internally to indicate one-arm or
    two-arm analysis.

  - `intent` character. Denotes current/historical status of treatment
    and control groups.

## Details

Historical borrowing is assessed separately within each arm by comparing
event probabilities. With `method = "fixed"`, estimate the probability q
that the current event probability is below the historical one. The
discount input is `p = 2 * min(q, 1 - q)` and the weight is
`alpha_max * W(p)` for the selected discount function W. This comparison
is not a posterior probability of equality or equivalence. Choosing
`fix_alpha = TRUE` prespecifies the weight instead.

Each supplied group requires scalar, finite integer event counts and
sample sizes, with \\N \geq 1\\ and \\0 \leq y \leq N\\. Both zero-event
and all-event groups are allowed. Omit an optional group by leaving both
inputs `NULL`; zero sample size does not denote an absent group. Inputs
and borrowing settings are checked before sampling, including when
weights are prescribed. Weibull parameters must be positive and finite
when used to estimate a weight.

In the second stage, posterior estimation is performed where the
discount function parameter, `alpha`, is incorporated in all posterior
estimation procedures.

To carry out a single arm (OPC) analysis, data for the current treatment
(`y_t` and `N_t`) and historical treatment (`y0_t` and `N0_t`) must be
input. The results are then based on the posterior distribution of the
current data augmented by the historical data.

To carry our a two-arm (RCT) analysis, data for the current treatment
and at least one of current or historical control data must be input.
The results are then based on the posterior distribution of the
difference between current treatment and control, augmented by available
historical data.

For more details, see the `bdpbinomial` vignette:  
[`vignette("bdpbinomial-vignette", package = "bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpbinomial-vignette.md)

## References

Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
Modeling and Simulation working group.(2017) Incorporation of stochastic
engineering models as prior information in Bayesian medical device
trials. *Journal of Biopharmaceutical Statistics*, 1-15.

## See also

[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpbinomial-method.md),
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpbinomial-method.md),
and
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpbinomial-method.md)
for details of each of the supported methods.

## Examples

``` r
# One-arm trial (OPC) example
fit <- bdpbinomial(
  y_t = 10,
  N_t = 500,
  y0_t = 25,
  N0_t = 250,
  method = "fixed"
)
summary(fit)
#> 
#>     One-armed bdp binomial
#> 
#> Current treatment data: 10 and 500
#> Historical treatment data: 25 and 250
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0
#> Discount function value (alpha) - treatment: 0
#> 95 percent equal-tailed credible interval: 
#>  0.0111  0.036
#> Posterior medians:
#>  0.0211
print(fit)
#> 
#>     One-armed bdp binomial
#> 
#> Current treatment data: 10 and 500
#> Historical treatment data: 25 and 250
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0
#> Discount function value (alpha) - treatment: 0
#> 95 percent equal-tailed credible interval: 
#>  0.0111  0.036
#> Posterior medians:
#>  0.0211
if (FALSE) { # \dontrun{
plot(fit)
} # }

# Two-arm (RCT) example
fit2 <- bdpbinomial(
  y_t = 10,
  N_t = 500,
  y0_t = 25,
  N0_t = 250,
  y_c = 8,
  N_c = 500,
  y0_c = 20,
  N0_c = 250,
  method = "fixed"
)
summary(fit2)
#> 
#>     Two-armed bdp binomial
#> 
#> Current treatment data: 10 and 500
#> Current control data: 8 and 500
#> Historical treatment data: 25 and 250
#> Historical control data: 20 and 250
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0
#> Discount function value (alpha) - treatment: 0
#> Stochastic comparison (p_hat) - control (current vs. historical data): 0
#> Discount function value (alpha) - control: 0
#> Treatment effect: treatment - control
#> 95 percent equal-tailed credible interval:
#>   -0.013 0.0215
#> Posterior medians:
#> treatment  control
#>   0.0214  0.0173
print(fit2)
#> 
#>     Two-armed bdp binomial
#> 
#> Current treatment data: 10 and 500
#> Current control data: 8 and 500
#> Historical treatment data: 25 and 250
#> Historical control data: 20 and 250
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0
#> Discount function value (alpha) - treatment: 0
#> Stochastic comparison (p_hat) - control (current vs. historical data): 0
#> Discount function value (alpha) - control: 0
#> Treatment effect: treatment - control
#> 95 percent equal-tailed credible interval:
#>   -0.013 0.0215
#> Posterior medians:
#> treatment  control
#>   0.0214  0.0173
if (FALSE) { # \dontrun{
plot(fit2)
} # }
```
