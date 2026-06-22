# Bayesian Discount Prior: Historical Data Weight (alpha)

`alpha_discount` can be used to estimate the weight applied to
historical data in the context of a one- or two-arm clinical trial.
`alpha_discount` is not used internally but is given for educational
purposes.

## Usage

``` r
alpha_discount(
  p_hat = NULL,
  discount_function = "weibull",
  alpha_max = 1,
  weibull_scale = 0.135,
  weibull_shape = 3
)
```

## Arguments

- p_hat:

  scalar. The posterior probability of a stochastic comparison. This
  value can be the output of `posterior_probability` or a value between
  0 and 1.

- discount_function:

  character. Specify the discount function to use. Currently supports
  `weibull`, `scaledweibull`, and `identity`. The discount function
  `scaledweibull` scales the output of the Weibull CDF to have a max
  value of 1. The `identity` discount function uses the posterior
  probability directly as the discount weight. Default value is
  "`weibull`".

- alpha_max:

  scalar. Maximum weight the discount function can apply. Default is 1.

- weibull_scale:

  scalar. Scale parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 0.135.

- weibull_shape:

  scalar. Shape parameter of the Weibull discount function used to
  compute alpha, the weight parameter of the historical data. Default
  value is 3.

## Value

`alpha_discount` returns an object of class "alpha_discount".

An object of class `alpha_discount` contains the following:

- `alpha_hat`:

  scalar. The historical data weight.

## Details

This function is not used internally but is given for educational
purposes. Given inputs `p_hat`, `discount_function`, `alpha_max`,
`weibull_shape`, and `weibull_scale` the output is the weight that would
be applied to historical data in the context of a one- or two-arm
clinical trial.

## References

Haddad, T., Himes, A., Thompson, L., Irony, T., Nair, R. MDIC Computer
Modeling and Simulation working group.(2017) Incorporation of stochastic
engineering models as prior information in Bayesian medical device
trials. *Journal of Biopharmaceutical Statistics*, 1-15.

## Examples

``` r
alpha_discount(0.5)
#> [1] 1

alpha_discount(0.5, discount_function = "identity")
#> [1] 0.5
```
