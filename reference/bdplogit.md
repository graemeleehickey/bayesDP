# Bayesian Discount Prior: Two-Arm Logistic Regression

`bdplogit` is used to estimate the treatment effect in the presence of
covariates using the logistic regression analysis implementation of the
Bayesian discount prior. `summary` and `print` methods are supported.
Currently, the function only supports a two-arm clinical trial where all
of current treatment, current control, historical treatment, and
historical control data are present

## Usage

``` r
bdplogit(
  formula = formula,
  data = data,
  data0 = NULL,
  prior_treatment_effect = NULL,
  prior_control_effect = NULL,
  prior_treatment_sd = NULL,
  prior_control_sd = NULL,
  prior_covariate_effect = 0,
  prior_covariate_sd = 10000,
  number_mcmc_alpha = 5000,
  number_mcmc_beta = 10000,
  discount_function = "identity",
  alpha_max = 1,
  fix_alpha = FALSE,
  weibull_scale = 0.135,
  weibull_shape = 3,
  method = "mc"
)
```

## Arguments

- formula:

  an object of class "formula." See "Details" for more information,
  including specification of treatment data indicators.

- data:

  a data frame containing the current data variables in the model. A
  column named `treatment` must be present; `treatment` must be binary
  and indicate treatment group vs. control group.

- data0:

  a data frame containing the historical data variables in the model.
  The column labels of data and data0 must match.

- prior_treatment_effect:

  scalar. The historical adjusted treatment effect. If left `NULL`,
  value is estimated from the historical data.

- prior_control_effect:

  scalar. The historical adjusted control effect. If left `NULL`, value
  is estimated from the historical data.

- prior_treatment_sd:

  scalar. The standard deviation of the historical adjusted treatment
  effect. If left `NULL`, value is estimated from the historical data.

- prior_control_sd:

  scalar. The standard deviation of the historical adjusted control
  effect. If left `NULL`, value is estimated from the historical data.

- prior_covariate_effect:

  vector. The prior mean(s) of the covariate effect(s). Default value is
  zero. If a single value is input, the the scalar is repeated to the
  length of the input covariates. Otherwise, care must be taken to
  ensure the length of the input matches the number of covariates.

- prior_covariate_sd:

  vector. The prior standard deviation(s) of the covariate effect(s).
  Default value is 1e4. If a single value is input, the the scalar is
  repeated to the length of the input covariates. Otherwise, care must
  be taken to ensure the length of the input matches the number of
  covariates. Note that covariate effects are intentionally assigned a
  near-zero historical discount weight, so their priors are effectively
  flat: the value supplied here has negligible influence on the
  posterior, and the effective prior standard deviation is much larger
  than the nominal value (by a factor of roughly 1e6 at the default). In
  practice the covariate coefficients are estimated almost entirely from
  the current data, and the posterior is insensitive to this argument
  across many orders of magnitude.

- number_mcmc_alpha:

  scalar. Number of Monte Carlo simulations for estimating the
  historical data weight. Default is 5000.

- number_mcmc_beta:

  scalar. Number of Monte Carlo simulations for estimating beta, the
  vector of regression coefficients. Default is 10000.

- discount_function:

  character. Specify the discount function to use. Currently supports
  `weibull`, `scaledweibull`, and `identity`. The discount function
  `scaledweibull` scales the output of the Weibull CDF to have a max
  value of 1. The `identity` discount function uses the posterior
  probability directly as the discount weight. Default value is
  "`identity`".

- alpha_max:

  scalar in \[0, 1\]. Maximum weight the discount function can apply.
  Default is 1. Users may specify a vector of two values where the first
  value is used to weight the historical treatment group and the second
  value is used to weight the historical control group.

- fix_alpha:

  logical. Fix alpha at alpha_max? Default value is FALSE. If TRUE,
  comparison models are not fitted and the prescribed weights are used
  directly. Historical-prior estimation and current-data model checks
  still apply.

- weibull_scale:

  scalar. Scale parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 0.135. Users may specify a vector of two values where the
  first value is used to estimate the weight of the historical treatment
  group and the second value is used to estimate the weight of the
  historical control group. Not used when `discount_function` =
  "identity".

- weibull_shape:

  scalar. Shape parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 3. Users may specify a vector of two values where the first
  value is used to estimate the weight of the historical treatment group
  and the second value is used to estimate the weight of the historical
  control group. Not used when `discount_function` = "identity".

- method:

  character. Set to `"fixed"`, the only supported method. It estimates
  each historical weight once and holds it fixed throughout the
  analysis, unless `fix_alpha = TRUE` prescribes the weights. The legacy
  default `"mc"` is not supported and raises an error.

## Value

`bdplogit` returns an object of class "bdplm", reusing its print and
summary methods. Residual summaries use observed binary outcomes minus
fitted probabilities from the posterior mean coefficients. The centering
means are saved in `args1$covariate_means`; older saved fits without
these means must be refitted before summarising residuals.

The returned object is a list containing at least the following
components:

- `posterior`:

  data frame. Posterior draws of the intercept, treatment effect, and
  covariate effects on the log-odds scale.

- `alpha_discount`:

  data frame. One row of estimated or prescribed historical weights,
  with columns `treatment` and `control`. These weights are not
  posterior probabilities of equality.

- `estimates`:

  list. Posterior means and posterior SDs of the intercept, treatment
  effect, and covariate effects.

## Details

`bdplogit` uses a two-stage approach for determining the strength of
historical data in estimation of an adjusted mean or covariate effect.
In the first stage, a *discount function* is used that defines the
maximum strength of the historical data and discounts based on
disagreement with the current data. Disagreement between current and
historical data is assessed by fitting a logistic regression separately
within each arm, with a historical indicator and the covariates, using
flat coefficient priors. If \\q\\ is the posterior probability that the
historical-indicator coefficient is positive, the comparison is \\p =
2\min(q, 1-q)\\. The discount function transforms \\p\\ into the
historical weight `alpha`. This is not a posterior probability of
equality or equivalence.

In the second stage, posterior estimation is performed where the
discount function parameter, `alpha`, is used to weight the historical
data effects.

The formula must include an intercept (i.e., do not use `-1` in the
formula) and both data and data0 must be present. The column names of
data and data0 must match. See `examples` below for example usage.

Formula offsets are not supported. Calls to
[`offset()`](https://rdrr.io/r/stats/offset.html) or
[`stats::offset()`](https://rdrr.io/r/stats/offset.html) cause an error
before fitting, including offsets that are zero or cancel each other.
Prescribing weights, disabling borrowing, or supplying historical priors
does not bypass this restriction.

Internally, the model is parameterized without a fixed intercept, using
separate treatment and control (log-odds) means rather than an intercept
plus a treatment effect. In this parameterization the arm means are the
fitted values at covariate = 0, so when covariates are not centered the
treatment- and control-mean estimators become strongly correlated and
their standard errors are extrapolation errors at covariate = 0. This
distorts the discount-prior construction and the historical borrowing.
To guard against this, `bdplogit` automatically mean-centers each
covariate on its pooled (current plus historical) mean before fitting,
and back-transforms the reported intercept onto the original covariate
scale. Specifically, if \\\bar{x}\_j\\ is the pooled mean of covariate
\\j\\ and \\\beta_j\\ its estimated effect, the reported intercept is
\\\beta_0 = \beta_0^{c} - \sum_j \beta_j \bar{x}\_j\\, where
\\\beta_0^{c}\\ is the control-arm log-odds estimated on the centered
scale. The treatment effect and covariate slopes are unchanged by
centering; the reported `intercept` remains the control-arm log-odds at
covariate = 0 on the original scale. As a result, estimates are
invariant to a location shift of any covariate, and users do not need to
center covariates themselves.

Posterior sampling uses
[`MCMCpack::MCMClogit()`](https://rdrr.io/pkg/MCMCpack/man/MCMClogit.html).
Before fitting, model matrices are checked for rank deficiency, and
complete or quasi-complete separation is detected using
`detectseparation`. These checks include separation caused by
covariates, as well as groups with no events or with all events. An
error identifies the affected stage:

- Each arm's comparison model must have a proper posterior under its
  flat coefficient priors. Separation or rank deficiency prevents
  estimation of that arm's borrowing weight.

- When any historical arm-prior mean or SD is estimated from the
  historical regression, that regression must have a full-rank design
  and no separation. Numerically large but finite estimates from a
  separated regression are not used as historical priors.

- The current-data regression must also have a full-rank design and no
  separation because the sampler uses an unpenalized logistic regression
  to construct starting values and proposals. This is a limitation of
  the supported sampling procedure, not a claim that separation makes a
  posterior with proper coefficient priors improper.

With `fix_alpha = TRUE`, the comparison models and their checks are
skipped; their internal `p_hat` values are `NA_real_`. Historical
regression checks are skipped only when all four arm-prior means and SDs
are supplied. Current-data checks always apply. Prescribing weights
alone therefore does not make separated datasets supported.

## Examples

``` r
# Set sample sizes
n_t <- 30 # Current treatment sample size
n_c <- 30 # Current control sample size
n_t0 <- 80 # Historical treatment sample size
n_c0 <- 80 # Historical control sample size

# Treatment group vectors for current and historical data
treatment <- c(rep(1, n_t), rep(0, n_c))
treatment0 <- c(rep(1, n_t0), rep(0, n_c0))

# Simulate a covariate effect for current and historical data
x <- rnorm(n_t + n_c, 1, 1)
x0 <- rnorm(n_t0 + n_c0, 1, 1)

# Simulate a binary outcome on the logit scale:
# - Intercept of -1 for current and historical data
# - Treatment (log-odds) effect of 1.0 for current data
# - Treatment (log-odds) effect of 0.9 for historical data
# - Covariate effect of 0.3 for current and historical data
Y <- rbinom(n_t + n_c, 1, plogis(-1 + 1.0 * treatment + 0.3 * x))
Y0 <- rbinom(n_t0 + n_c0, 1, plogis(-1 + 0.9 * treatment0 + 0.3 * x0))

# Place data into separate treatment and control data frames and
# assign historical = 0 (current) or historical = 1 (historical)
df_ <- data.frame(Y = Y, treatment = treatment, x = x)
df0 <- data.frame(Y = Y0, treatment = treatment0, x = x0)

# Fit model using default settings
fit <- bdplogit(
  formula = Y ~ treatment + x, data = df_, data0 = df0,
  method = "fixed"
)

# Look at estimates and discount weight
summary(fit)
#> 
#> Call:
#> bdplogit(formula = Y ~ treatment + x, data = df_, data0 = df0, 
#>     method = "fixed")
#> 
#> Residuals:
#>     Min     1Q Median    3Q   Max
#>  -0.591 -0.378 -0.265 0.491 0.723
#> 
#> Coefficients:
#>             Estimate Std. Error
#> (Intercept)  -1.0822     0.3901
#> treatment     0.8140     0.4008
#> x             0.2548     0.2830
#> 
#> Discount function value (alpha):
#>  treatment control
#>      0.244   0.422
#> 
print(fit)
#> 
#> Call:
#> bdplogit(formula = Y ~ treatment + x, data = df_, data0 = df0, 
#>     method = "fixed")
#> 
#> 
#> Coefficients:
#>  (Intercept) treatment     x
#>       -1.082     0.814 0.255
#> 
#> 
#> Discount function value (alpha):
#>  treatment control
#>      0.244   0.422
#> 
```
