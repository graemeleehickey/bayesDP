# Bayesian Discount Prior: Two-Arm Linear Regression

`bdplm` is used to estimate the treatment effect in the presence of
covariates using the regression analysis implementation of the Bayesian
discount prior. `summary` and `print` methods are supported. Currently,
the function only supports a two-arm clinical trial where all of current
treatment, current control, historical treatment, and historical control
data are present

## Usage

``` r
bdplm(
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
  number_mcmc_sigmagrid = 5000,
  number_mcmc_sigma = 100,
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

- number_mcmc_sigmagrid:

  scalar. Grid size for computing sigma. Default is 5000. See "Details"
  for more information.

- number_mcmc_sigma:

  scalar. Number of Monte Carlo simulations for estimating sigma.
  Default is 1000. See "Details" for more information.

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

  scalar. Maximum weight the discount function can apply. Default is 1.
  Users may specify a vector of two values where the first value is used
  to weight the historical treatment group and the second value is used
  to weight the historical control group.

- fix_alpha:

  logical. Fix alpha at alpha_max? Default value is FALSE.

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

  character. Analysis method with respect to estimation of the weight
  paramter alpha. Default method "`mc`" estimates alpha for each Monte
  Carlo iteration. Alternate value "`fixed`" estimates alpha once and
  holds it fixed throughout the analysis. See the the `bdplm` vignette  
  [`vignette("bdplm-vignette", package="bayesDP")`](https://graemeleehickey.github.io/bayesDP/articles/bdplm-vignette.md)
  for more details.

## Value

`bdplm` returns an object of class "bdplm".

An object of class "`bdplm`" is a list containing at least the following
components:

- `posterior`:

  data frame. The posterior draws of the covariates, the intercept, and
  the treatment effect. The grid of sigma values are included.

- `alpha_discount`:

  vector. The posterior probability of the stochastic comparison between
  the current and historical data for each of the treatment and control
  arms. If `method="mc"`, the result is a matrix of estimates, otherwise
  for `method="fixed"`, the result is a vector of estimates.

- `estimates`:

  list. The posterior means and standard errors of the intercept,
  treatment effect, covariate effect(s) and error variance.

## Details

`bdplm` uses a two-stage approach for determining the strength of
historical data in estimation of an adjusted mean or covariate effect.
In the first stage, a *discount function* is used that that defines the
maximum strength of the historical data and discounts based on
disagreement with the current data. Disagreement between current and
historical data is determined by stochastically comparing the respective
posterior distributions under noninformative priors. Here with a two-arm
regression analysis, the comparison is the proability (`p`) that the
covariate effect of an historical data indicator is significantly
different from zero. The comparison metric `p` is then input into the
discount function and the final strength of the historical data is
returned (`alpha`).

In the second stage, posterior estimation is performed where the
discount function parameter, `alpha`, is used to weight the historical
data effects.

The formula must include an intercept (i.e., do not use `-1` in the
formula) and both data and data0 must be present. The column names of
data and data0 must match. See `examples` below for example usage.

Internally, the model is parameterized without a fixed intercept, using
separate treatment and control means rather than an intercept plus a
treatment effect. In this parameterization the arm means are the fitted
values at covariate = 0, so when covariates are not centered the
treatment- and control-mean estimators become strongly correlated and
their standard errors are extrapolation errors at covariate = 0. This
distorts the discount-prior construction and the historical borrowing.
To guard against this, `bdplm` automatically mean-centers each covariate
on its pooled (current plus historical) mean before fitting, and
back-transforms the reported intercept onto the original covariate
scale. Specifically, if \\\bar{x}\_j\\ is the pooled mean of covariate
\\j\\ and \\\beta_j\\ its estimated effect, the reported intercept is
\\\beta_0 = \beta_0^{c} - \sum_j \beta_j \bar{x}\_j\\, where
\\\beta_0^{c}\\ is the control-arm mean estimated on the centered scale.
The treatment effect and covariate slopes are unchanged by centering;
the reported `intercept` remains the control-arm mean at covariate = 0
on the original scale. As a result, estimates are invariant to a
location shift of any covariate, and users do not need to center
covariates themselves.

The underlying model results in a marginal posterior distribution for
the error variance `sigma2` that does not have a known distribution.
Thus, we use a grid search to draw samples of `sigma2`. First, the
marginal posterior is evaluated at a grid of `number_mcmc_sigmagrid`
`sigma2` values. The bounds of the grid are taken as approximately six
standard deviations from an estimate of `sigma2` using the `lm`
function. Next, `number_mcmc_sigma` posterior draws of `sigma2` are made
by sampling with replacement from the grid with each value having the
corresponding marginal posterior probability of being selected. With
posterior draws of `sigma2` in hand, we can make posterior draws of the
regression coefficients.

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
x <- rnorm(n_t + n_c, 1, 5)
x0 <- rnorm(n_t0 + n_c0, 1, 5)

# Simulate outcome:
# - Intercept of 10 for current and historical data
# - Treatment effect of 31 for current data
# - Treatment effect of 30 for historical data
# - Covariate effect of 3 for current and historical data
Y <- 10 + 31 * treatment + x * 3 + rnorm(n_t + n_c, 0, 5)
Y0 <- 10 + 30 * treatment0 + x0 * 3 + rnorm(n_t0 + n_c0, 0, 5)

# Place data into separate treatment and control data frames and
# assign historical = 0 (current) or historical = 1 (historical)
df_ <- data.frame(Y = Y, treatment = treatment, x = x)
df0 <- data.frame(Y = Y0, treatment = treatment0, x = x0)

# Fit model using default settings
fit <- bdplm(
  formula = Y ~ treatment + x, data = df_, data0 = df0,
  method = "fixed"
)

# Look at estimates and discount weight
summary(fit)
#> 
#> Call:
#> bdplm(formula = Y ~ treatment + x, data = df_, data0 = df0, method = "fixed")
#> 
#> Residuals:
#>     Min     1Q Median    3Q    Max
#>  -7.096 -0.492  2.569 7.794 15.307
#> 
#> Coefficients:
#>             Estimate Std. Error
#> (Intercept)   9.7721     0.7334
#> treatment    30.3256     1.0342
#> x             2.8382     0.1455
#> 
#> Discount function value (alpha):
#>  treatment control
#>     0.2776  0.3348
#> 
#> Residual standard error: 5.681
print(fit)
#> 
#> Call:
#> bdplm(formula = Y ~ treatment + x, data = df_, data0 = df0, method = "fixed")
#> 
#> 
#> Coefficients:
#>  (Intercept) treatment     x
#>        9.772    30.326 2.838
#> 
#> 
#> Discount function value (alpha):
#>  treatment control
#>     0.2776  0.3348
#> 
```
