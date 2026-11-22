# Bayesian Discount Prior: Survival Analysis

`bdpsurvival` is used to estimate the survival probability (single arm
trial; OPC) or hazard ratio (two-arm trial; RCT) for right-censored data
using the survival analysis implementation of the Bayesian discount
prior. Current treatment data are required; control data in a two-arm
analysis can be current, historical, or both. This code is modeled after
the methodologies developed in Haddad et al. (2017).

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
  for more information. Response expressions such as
  `Surv(followup / 12, status == 1)` and
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) are
  supported. Only right-censored responses are supported.

- data:

  a data frame containing the current data variables in the model.
  Supply the variables used in the formula. Evaluated times must be
  finite and non-negative, with no missing outcomes. If treatment is
  specified, use numeric 0/1 or logical values, with no missing values;
  factors and character codes are not converted automatically.

- data0:

  optional. A data frame containing the historical data variables in the
  model, subject to the same response and treatment checks as current
  data. It must supply the variables used in the same formula.

- breaks:

  numeric vector. Positive, finite internal cutpoints of the piecewise
  exponential model, without zero or infinity. Values are sorted and
  duplicates removed. Use `numeric(0)` for one constant hazard. Default
  cutpoints are pooled observed-time quantiles.

- a0:

  positive, finite scalar. Prior value for the gamma shape of the
  piecewise exponential hazards. Default is 0.1.

- b0:

  positive, finite scalar. Prior value for the gamma rate of the
  piecewise exponential hazards. Default is 0.1.

- surv_time:

  positive, finite scalar. Time of interest for computing the
  probability of survival for a single arm (OPC) trial. Default is
  overall, i.e., current+historical, median observed event-or-censoring
  time.

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

- number_mcmc:

  integer of at least 2. Number of Monte Carlo simulations. Default is
  10000.

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
  Carlo iteration, giving a distribution of generated weights.
  Increasing the number of draws improves numerical precision but does
  not remove the spread of that distribution. These weights are not
  draws from a separately specified prior/posterior model for alpha.
  Alternate value "`fixed`" estimates alpha once and holds it fixed
  throughout the analysis. See the `bdpsurvival` vignette  
  [`vignette("bdpsurvival-vignette", package = "bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpsurvival-vignette.md)
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

  - `p_hat` numeric. The two-sided current-historical comparison in \[0,
    1\]: one value for `method = "fixed"`, or `number_mcmc` values for
    `method = "mc"`. See Details for the transformation.

  - `posterior_survival` matrix. If one-arm trial, one column with
    `number_mcmc` rows containing the posterior probability draws of
    survival at `surv_time`.

  - `posterior_flat_survival` matrix. If one-arm trial, one column with
    `number_mcmc` rows containing the probability draws of survival at
    `surv_time` for the current treatment not augmented by historical
    treatment.

  - `prior_survival` matrix. If one-arm trial, one column with
    `number_mcmc` rows containing the probability draws of survival at
    `surv_time` for the historical treatment.

  - `posterior_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks) + 1` columns containing the posterior draws of the
    piecewise hazards for each interval break point.

  - `posterior_flat_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks) + 1` columns containing the draws of piecewise
    hazards for each interval break point for current treatment not
    augmented by historical treatment.

  - `prior_hazard` matrix. A matrix with `number_mcmc` rows and
    `length(breaks) + 1` columns containing the draws of piecewise
    hazards for each interval break point for historical treatment.

  - `posterior_log_hazard`, `posterior_flat_log_hazard`,
    `prior_log_hazard` matrices. Log hazards corresponding to the three
    hazard matrices. These retain finite log draws when very small
    hazards round to zero and are used for borrowing comparisons and the
    two-arm effect.

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

  - `participants` data frame. One record per participant, with observed
    `time`, event `status`, `treatment`, and `historical` source. Used
    for reported sample sizes and risk sets; `data` contains the split
    interval records used for hazard estimation.

## Details

Historical borrowing is assessed separately within each arm by comparing
survival outcomes. With `method = "fixed"`, estimate a directional
probability q comparing current and historical outcomes within an arm.
The discount input is `p = 2 * min(q, 1 - q)` and the weight is
`alpha_max * W(p)` for the selected discount function W. This comparison
is not a posterior probability of equality or equivalence. Choosing
`fix_alpha = TRUE` prespecifies the weight instead.

For one arm, q compares survival at `surv_time`; for two arms it
compares current and historical hazards separately within treatment and
control using analytic inverse-variance weights for the interval log
hazard differences. With `method = "mc"`, each draw's difference is
divided by its posterior standard deviation and mapped to
`2 * pnorm(abs(Z), lower.tail = FALSE)`. The variance is calculated from
exact Gamma moments: hazard variances for log survival differences, and
trigamma log-hazard variances for two-arm comparisons. The normal-tail
mapping remains an approximation, especially with sparse events. The
returned `p_hat` uses the selected method's transformation even with
`fix_alpha = TRUE`; the prescribed weight is unaffected.

In the second stage, the historical likelihood is raised to the weight
`alpha`, and the baseline Gamma prior is included once. For an interval
with current event count D and exposure T, and historical event count D0
and exposure T0, the augmented Gamma shape and rate are
`a0 + D + alpha * D0` and `b0 + T + alpha * T0`. The prior shape `a0` is
distinct from the historical weight alpha. At alpha = 0 this gives the
current-only posterior; at alpha = 1 it gives the pooled posterior using
the same interval boundaries and baseline prior. With `method = "mc"`,
each draw uses its generated weight; the resulting marginal distribution
is a mixture of Gamma distributions.

To carry out a single arm (OPC) analysis, data for the current and
historical treatments are specified in separate data frames, data and
data0, respectively. The formula is evaluated separately in each data
frame to obtain the survival (censor) time and event indicator. Missing
outcomes or treatment values cause an error rather than dropping rows.
Counting-process, interval-censored, and multistate responses are not
supported. The results are then based on the posterior probability of
survival at `surv_time` for the current data augmented by the historical
data.

Two-arm (RCT) analyses use the same response specification, with an
additional column named 'treatment' to denote treatment and control
data. Use 1 (or TRUE) for treatment and 0 (or FALSE) for control in each
data frame. The current data are augmented by historical data (if
present) and the results are then based on the posterior distribution of
the hazard ratio between the treatment and control groups.

For more details, see the `bdpsurvival` vignette:  
[`vignette("bdpsurvival-vignette", package = "bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdpsurvival-vignette.md)

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
#>   n events surv_time median lower 95% CrI upper 95% CrI
#>  50     31         5 0.7756        0.6608        0.8695
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
#>   Current treatment: n = 50, number of events = 35
#>   Current control: n = 50, number of events = 36
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.1478
#> Stochastic comparison (p_hat) - control (current vs. historical data): 4e-04
#> Discount function value (alpha) - treatment: 0.1478
#> Discount function value (alpha) - control: 4e-04
#> 
#>              coef exp(coef) se(coef) lower 95% CI upper 95% CI
#> treatment -0.5685    0.5664   0.2453      -1.0563      -0.0889

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
#>   Current treatment: n = 50, number of events = 35
#>   Current control: n = 50, number of events = 36
#> Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.1592
#> Stochastic comparison (p_hat) - control (current vs. historical data): 0.0012
#> Discount function value (alpha) - treatment: 1
#> Discount function value (alpha) - control: 1
#> 
#>              coef exp(coef) se(coef) lower 95% CI upper 95% CI
#> treatment -1.0445    0.3519   0.1837      -1.4078      -0.6915
```
