# Bayesian Discount Prior: Survival Analysis

`bdpsurvival` is used to estimate the survival probability (single arm
trial; OPC) or hazard ratio (two-arm trial; RCT) for right-censored data
using the survival analysis implementation of the Bayesian discount
prior. In the current implementation, a two-arm analysis requires all of
current treatment, current control, historical treatment, and historical
control data. This code is modeled after the methodologies developed in
Haddad et al. (2017).

## Usage

``` r
bdpsurvival(
  formula = formula,
  data = data,
  data0 = NULL,
  breaks = NULL,
  a0 = 0.1,
  b0 = 0.1,
  surv_time = NULL,
  discount_function = "identity",
  alpha_max = 1,
  fix_alpha = FALSE,
  number_mcmc = 10000,
  weibull_scale = 0.135,
  weibull_shape = 3,
  method = "mc",
  compare = TRUE
)
```

## Arguments

- formula:

  an object of class "formula." Must have a survival object on the left
  side and at most one input on the right side, treatment. See "Details"
  for more information.

- data:

  a data frame containing the current data variables in the model.
  Columns denoting 'time' and 'status' must be present. See "Details"
  for required structure.

- data0:

  optional. A data frame containing the historical data variables in the
  model. If present, the column labels of data and data0 must match.

- breaks:

  vector. Breaks (interval starts) used to compose the breaks of the
  piecewise exponential model. Do not include zero. Default breaks are
  the quantiles of the input times.

- a0:

  scalar. Prior value for the gamma shape of the piecewise exponential
  hazards. Default is 0.1.

- b0:

  scalar. Prior value for the gamma rate of the piecewise exponential
  hazards. Default is 0.1.

- surv_time:

  scalar. Survival time of interest for computing the probability of
  survival for a single arm (OPC) trial. Default is overall, i.e.,
  current+historical, median survival time.

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

- number_mcmc:

  scalar. Number of Monte Carlo simulations. Default is 10000.

- weibull_scale:

  scalar. Scale parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 0.135. For a two-arm trial, users may specify a vector of two
  values where the first value is used to estimate the weight of the
  historical treatment group and the second value is used to estimate
  the weight of the historical control group.

- weibull_shape:

  scalar. Shape parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 3. For a two-arm trial, users may specify a vector of two
  values where the first value is used to estimate the weight of the
  historical treatment group and the second value is used to estimate
  the weight of the historical control group.

- method:

  character. Analysis method with respect to estimation of the weight
  parameter alpha. Default method "`mc`" estimates alpha for each Monte
  Carlo iteration. Alternate value "`fixed`" estimates alpha once and
  holds it fixed throughout the analysis. See the the `bdpsurvival`
  vignette  
  [`vignette("bdpsurvival-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpsurvival-vignette.md)
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

`bdpsurvival` returns an object of class "bdpsurvival". The functions
[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpsurvival-method.md)
and
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpsurvival-method.md)
are used to obtain and print a summary of the results, including user
inputs. The
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpsurvival-method.md)
function displays visual outputs as well.

An object of class "`bdpsurvival`" is a list containing at least the
following components:

- `posterior_treatment`:

  list. Entries contain values related to the treatment group:

  - `alpha_discount` numeric. Alpha value, the weighting parameter of
    the historical data.

  - `p_hat` numeric. The posterior probability of the stochastic
    comparison between the current and historical data.

  - `posterior_survival` vector. If one-arm trial, a vector of length
    `number_mcmc` containing the posterior probability draws of survival
    at `surv_time`.

  - `posterior_flat_survival` vector. If one-arm trial, a vector of
    length `number_mcmc` containing the probability draws of survival at
    `surv_time` for the current treatment not augmented by historical
    treatment.

  - `prior_survival` vector. If one-arm trial, a vector of length
    `number_mcmc` containing the probability draws of survival at
    `surv_time` for the historical treatment.

  - `posterior_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks)` columns containing the posterior draws of the
    piecewise hazards for each interval break point.

  - `posterior_flat_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks)` columns containing the draws of piecewise hazards
    for each interval break point for current treatment not augmented by
    historical treatment.

  - `prior_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks)` columns containing the draws of piecewise hazards
    for each interval break point for historical treatment.

- `posterior_control`:

  list. If two-arm trial, contains values related to the control group
  analogous to the `posterior_treatment` output.

- `final`:

  list. Contains the final comparison object, dependent on the analysis
  type:

  - One-arm analysis: vector. Posterior chain of survival probability at
    requested time.

  - Two-arm analysis: vector. Posterior chain of log-hazard rate
    comparing treatment and control groups.

- `args1`:

  list. Entries contain user inputs. In addition, the following elements
  are output:

  - `S_t`, `S_c`, `S0_t`, `S0_c` survival objects. Used internally to
    pass survival data between functions.

  - `arm2` logical. Used internally to indicate one-arm or two-arm
    analysis.

## Details

`bdpsurvival` uses a two-stage approach for determining the strength of
historical data in estimation of a survival probability outcome. In the
first stage, a *discount function* is used that that defines the maximum
strength of the historical data and discounts based on disagreement with
the current data. Disagreement between current and historical data is
determined by stochastically comparing the respective posterior
distributions under noninformative priors. With a single arm survival
data analysis, the comparison is the probability (`p`) that the current
survival is less than the historical survival. For a two-arm survival
data, analysis the comparison is the probability that the hazard ratio
comparing treatment and control is different from zero. The comparison
metric `p` is then input into the discount function and the final
strength of the historical data is returned (alpha).

In the second stage, posterior estimation is performed where the
discount function parameter, `alpha`, is used incorporated in all
posterior estimation procedures.

To carry out a single arm (OPC) analysis, data for the current and
historical treatments are specified in separate data frames, data and
data0, respectively. The data frames must have matching columns denoting
time and status. The 'time' column is the survival (censor) time of the
event and the 'status' column is the event indicator. The results are
then based on the posterior probability of survival at `surv_time` for
the current data augmented by the historical data.

Two-arm (RCT) analyses are specified similarly to a single arm trial.
Again the input data frames must have columns denoting time and status,
but now an additional column named 'treatment' is required to denote
treatment and control data. The 'treatment' column must use 0 to
indicate the control group. The current data are augmented by historical
data (if present) and the results are then based on the posterior
distribution of the hazard ratio between the treatment and control
groups.

For more details, see the `bdpsurvival` vignette:  
[`vignette("bdpsurvival-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpsurvival-vignette.md)

## References

Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
Modeling and Simulation working group.(2017) Incorporation of stochastic
engineering models as prior information in Bayesian medical device
trials. *Journal of Biopharmaceutical Statistics*, 1-15.

## See also

[`summary`](https://graemeleehickey.github.io/bayesDP/reference/summary-bdpsurvival-method.md),
[`print`](https://graemeleehickey.github.io/bayesDP/reference/print-bdpsurvival-method.md),
and
[`plot`](https://graemeleehickey.github.io/bayesDP/reference/plot-bdpsurvival-method.md)
for details of each of the supported methods.

## Examples

``` r
# One-arm trial (OPC) example - survival probability at 5 years

# Collect data into data frames
df_ <- data.frame(
  status = rexp(50, rate = 1 / 30),
  time = rexp(50, rate = 1 / 20)
)
df_$status <- ifelse(df_$time < df_$status, 1, 0)

df0 <- data.frame(
  status = rexp(50, rate = 1 / 30),
  time = rexp(50, rate = 1 / 10)
)
df0$status <- ifelse(df0$time < df0$status, 1, 0)


fit1 <- bdpsurvival(Surv(time, status) ~ 1,
  data = df_,
  data0 = df0,
  surv_time = 5,
  method = "fixed"
)

print(fit1)
#> 
#>     One-armed bdp survival
#> 
#> 
#>   n events surv_time median lower 95% CI upper 95% CI
#>  50     29         5  0.794       0.6771       0.8832
if (FALSE) { # \dontrun{
plot(fit1)
} # }

# Two-arm trial example
# Collect data into data frames
df_ <- data.frame(
  time = c(
    rexp(50, rate = 1 / 20), # Current treatment
    rexp(50, rate = 1 / 10)
  ), # Current control
  status = rexp(100, rate = 1 / 40),
  treatment = c(rep(1, 50), rep(0, 50))
)
df_$status <- ifelse(df_$time < df_$status, 1, 0)

df0 <- data.frame(
  time = c(
    rexp(50, rate = 1 / 30), # Historical treatment
    rexp(50, rate = 1 / 5)
  ), # Historical control
  status = rexp(100, rate = 1 / 40),
  treatment = c(rep(1, 50), rep(0, 50))
)
df0$status <- ifelse(df0$time < df0$status, 1, 0)

fit2 <- bdpsurvival(Surv(time, status) ~ treatment,
  data = df_,
  data0 = df0,
  method = "fixed"
)
summary(fit2)
#> 
#>     Two-armed bdp survival
#> 
#> data:
#>   Current treatment: n = 181, number of events = 29
#>   Current control: n = 132, number of events = 40
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.633
#> Stochastic comparison (p_hat) - control (current vs. historical data): 0.0018
#> Discount function value (alpha) - treatment: 0.633
#> Discount function value (alpha) - control: 0.0018
#> 
#>              coef exp(coef) se(coef) lower 95% CI upper 95% CI
#> treatment -1.0753    0.3412   0.2337      -1.5287      -0.6158

### Fix alpha at 1
fit2_1 <- bdpsurvival(Surv(time, status) ~ treatment,
  data = df_,
  data0 = df0,
  fix_alpha = TRUE,
  method = "fixed"
)
summary(fit2_1)
#> 
#>     Two-armed bdp survival
#> 
#> data:
#>   Current treatment: n = 181, number of events = 29
#>   Current control: n = 132, number of events = 40
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.6106
#> Stochastic comparison (p_hat) - control (current vs. historical data): 8e-04
#> Discount function value (alpha) - treatment: 1
#> Discount function value (alpha) - control: 1
#> 
#>              coef exp(coef) se(coef) lower 95% CI upper 95% CI
#> treatment -1.3787    0.2519   0.1895      -1.7559      -1.0089
```
