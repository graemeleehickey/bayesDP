# Survival outcomes with bayesDP

## Purpose and scope

[`bdpsurvival()`](https://graemeleehickey.github.io/bayesDP/reference/bdpsurvival.md)
uses piecewise constant hazards for right-censored outcomes. In one arm,
the estimand is survival probability at a specified time. With two arms,
the reported effect is an approximate common log hazard ratio, treatment
relative to control.

Within each arm, historical events and exposure contribute through a
borrowing weight, while the baseline Gamma prior is included once. These
examples illustrate the calculations and their sensitivity to the
borrowing method; they do not establish operating characteristics for a
trial design.

See the [borrowing
overview](https://graemeleehickey.github.io/bayesDP/articles/borrowing.md)
for weighting choices and interpretation. The main examples use
`method = "fixed"` for exposition and then compare it with the default
`method = "mc"`.

## Data and follow-up

Supply current and historical data frames with compatible columns. Use
`Surv(time, status) ~ 1` for a one-arm analysis, or
`Surv(time, status) ~ treatment` for two arms. Code `status = 1` for an
event and `status = 0` for right censoring; `time` is the observed event
or censoring time. For two arms, code `treatment = 1` and control as 0,
and provide both arms in both studies for the examples below. Additional
covariate adjustment is not implemented.

Choose a time unit, clinically relevant survival horizon, and interval
boundaries. `breaks` contains positive internal cutpoints, without zero
or infinity. The final interval extends beyond the last cutpoint;
`breaks = numeric(0)` specifies one constant hazard. Cutpoints are
sorted and duplicates removed.

By default the cutpoints are the pooled observed-time 20th, 40th, 60th,
and 80th percentiles; these are quantiles of event-or-censoring times,
not estimated survival quantiles. The examples prespecify cutpoints and
the horizon.

## One-arm example with censoring

We simulate event times with mean 10 months for 60 current participants
and mean 11 months for 100 historical participants. Independent
censoring times have mean 18 months, with additional administrative
censoring at 18 months. The observed time is the minimum of event and
censoring times.

``` r

simulate_followup <- function(n, mean_event, mean_censor, horizon) {
  event_time <- rexp(n, rate = 1 / mean_event)
  censor_time <- pmin(rexp(n, rate = 1 / mean_censor), horizon)
  data.frame(time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time))
}
set.seed(42)
current <- simulate_followup(60, mean_event = 10, mean_censor = 18, horizon = 18)
historical <- simulate_followup(100, mean_event = 11, mean_censor = 18, horizon = 18)
knitr::kable(data.frame(
  source = c("Current", "Historical"),
  participants = c(nrow(current), nrow(historical)),
  events = c(sum(current$status), sum(historical$status)),
  censored = c(sum(current$status == 0), sum(historical$status == 0))
))
```

| source     | participants | events | censored |
|:-----------|-------------:|-------:|---------:|
| Current    |           60 |     36 |       24 |
| Historical |          100 |     54 |       46 |

``` r

set.seed(42)
current_only <- bdpsurvival(Surv(time, status) ~ 1,
  data = current, data0 = historical, breaks = c(3, 6, 12), surv_time = 6,
  a0 = 0.1, b0 = 0.1, method = "fixed", fix_alpha = TRUE, alpha_max = 0)
set.seed(42)
discounted <- bdpsurvival(Surv(time, status) ~ 1,
  data = current, data0 = historical, breaks = c(3, 6, 12), surv_time = 6,
  a0 = 0.1, b0 = 0.1, method = "fixed", fix_alpha = FALSE, alpha_max = 1)
set.seed(42)
full_borrowing <- bdpsurvival(Surv(time, status) ~ 1,
  data = current, data0 = historical, breaks = c(3, 6, 12), surv_time = 6,
  a0 = 0.1, b0 = 0.1, method = "fixed", fix_alpha = TRUE, alpha_max = 1)
fits <- list("Current only" = current_only, "Discounted" = discounted,
  "Full borrowing" = full_borrowing)
comparison <- do.call(rbind, lapply(fits, function(fit) {
  arm <- fit$posterior_treatment
  c(alpha = arm$alpha_discount, median = median(arm$posterior_survival),
    lower = unname(quantile(arm$posterior_survival, 0.025)),
    upper = unname(quantile(arm$posterior_survival, 0.975)))
}))
knitr::kable(comparison, digits = 4,
  col.names = c("Historical weight", "Median S(6)", "Lower 95% CrI", "Upper 95% CrI"))
```

|                | Historical weight | Median S(6) | Lower 95% CrI | Upper 95% CrI |
|:---------------|------------------:|------------:|--------------:|--------------:|
| Current only   |            0.0000 |      0.4945 |        0.3572 |        0.6294 |
| Discounted     |            0.0272 |      0.5028 |        0.3690 |        0.6366 |
| Full borrowing |            1.0000 |      0.6175 |        0.5342 |        0.6969 |

The estimated historical weight is 0.0272. The augmented posterior
median six-month survival is 0.5028. These summarize survival
probability at six months, not median survival time.

The current-only result uses the current participants and the specified
baseline prior. Full borrowing gives the posterior from pooled current
and historical participants under the same piecewise exponential model
and baseline prior. Keep the interval boundaries and survival horizon
identical when making these comparisons; default cutpoints depend on the
supplied data.

[`print()`](https://rdrr.io/r/base/print.html) gives a compact summary
at the requested horizon. Its sample size and event count refer to
current participants. The longer
[`summary()`](https://rdrr.io/r/base/summary.html) table reports risk
sets immediately before each observed event or censoring time, alongside
posterior median survival, posterior SDs, and 95% equal-tailed credible
intervals. Splitting follow-up into hazard intervals does not create
additional participants or observation times. These risk sets describe
the observed current cohort; historical borrowing changes the survival
estimates, not the numbers at risk.

``` r

print(discounted)
```

    ## 
    ##     One-armed bdp survival
    ## 
    ## 
    ##   n events surv_time median lower 95% CrI upper 95% CrI
    ##  60     36         6 0.5028         0.369        0.6366

``` r

plot(discounted, type = "survival", print = FALSE) +
  labs(x = "Follow-up (months)")
```

![Posterior median survival curves for the estimated-borrowing fit.
These are model-based curves, not Kaplan–Meier estimates; uncertainty is
summarized in the six-month
table.](bdpsurvival-vignette_files/figure-html/survival-curves-1.png)

Posterior median survival curves for the estimated-borrowing fit. These
are model-based curves, not Kaplan–Meier estimates; uncertainty is
summarized in the six-month table.

### Sensitivity to the borrowing method

Keep the data, priors, intervals, horizon, discount function, and
maximum weight the same, and change the borrowing procedure. For the
draw-specific method, summarize the distribution of weights as well as
the survival estimate. The weight percentiles describe variation
introduced by this procedure; they are not a credible interval from a
separate posterior model for the weight.

``` r

set.seed(42)
discounted_mc <- bdpsurvival(Surv(time, status) ~ 1,
  data = current, data0 = historical, breaks = c(3, 6, 12), surv_time = 6,
  a0 = 0.1, b0 = 0.1, method = "mc", fix_alpha = FALSE,
  discount_function = "identity", alpha_max = 1)
format_interval <- function(x) paste(sprintf("%.3f", quantile(x, c(.025, .975))), collapse = " to ")
methods_one <- list("Fixed" = discounted, "Draw-specific" = discounted_mc)
knitr::kable(do.call(rbind, lapply(names(methods_one), function(method) {
  arm <- methods_one[[method]]$posterior_treatment
  data.frame(Method = method, `Median weight` = median(arm$alpha_discount),
    `Weight 2.5–97.5%` = format_interval(arm$alpha_discount),
    `Median S(6)` = median(arm$posterior_survival),
    `95% CrI for S(6)` = format_interval(arm$posterior_survival), check.names = FALSE)
})), digits = 3)
```

| Method        | Median weight | Weight 2.5–97.5% | Median S(6) | 95% CrI for S(6) |
|:--------------|--------------:|:-----------------|------------:|:-----------------|
| Fixed         |         0.027 | 0.027 to 0.027   |       0.503 | 0.369 to 0.637   |
| Draw-specific |         0.040 | 0.000 to 0.783   |       0.521 | 0.380 to 0.658   |

The two procedures need not give the same mean or median weight, even
with many simulation draws. Their definitions differ: one transforms an
estimated directional probability, whereas the other transforms each
standardized difference. Neither weight is the posterior probability of
equivalence between studies.

## Two-arm example

Add independent current and historical controls with mean event time 12
months, using the same censoring mechanism. All four sources are
present.

``` r

set.seed(101)
control <- simulate_followup(60, mean_event = 12, mean_censor = 18, horizon = 18)
control0 <- simulate_followup(100, mean_event = 12, mean_censor = 18, horizon = 18)
current_two <- rbind(transform(current, treatment = 1), transform(control, treatment = 0))
historical_two <- rbind(transform(historical, treatment = 1), transform(control0, treatment = 0))
set.seed(42)
fit_two <- bdpsurvival(Surv(time, status) ~ treatment,
  data = current_two, data0 = historical_two, breaks = c(3, 6, 12),
  a0 = 0.1, b0 = 0.1, method = "fixed", fix_alpha = FALSE, alpha_max = 1)
summary(fit_two)
```

    ## 
    ##     Two-armed bdp survival
    ## 
    ## data:
    ##   Current treatment: n = 60, number of events = 36
    ##   Current control: n = 60, number of events = 30
    ## Stochastic comparison (p_hat) - treatment (current vs. historical data): 0.1116
    ## Stochastic comparison (p_hat) - control (current vs. historical data): 0.8042
    ## Discount function value (alpha) - treatment: 0.1116
    ## Discount function value (alpha) - control: 0.8042
    ## 
    ##             coef exp(coef) se(coef) lower 95% CI upper 95% CI
    ## treatment 0.0621    1.0641   0.1955      -0.3275       0.4367

``` r

log_hr <- as.numeric(fit_two$final$posterior_loghazard)
knitr::kable(rbind(
  "Log hazard ratio" = quantile(log_hr, c(0.025, 0.5, 0.975)),
  "Hazard ratio" = quantile(exp(log_hr), c(0.025, 0.5, 0.975))
), digits = 4)
```

|                  |    2.5% |    50% |  97.5% |
|:-----------------|--------:|-------:|-------:|
| Log hazard ratio | -0.3275 | 0.0653 | 0.4367 |
| Hazard ratio     |  0.7208 | 1.0675 | 1.5476 |

The summary counts 60 current treatment participants and 60 current
control participants, regardless of how many intervals their follow-up
spans. The treatment and control historical weights are 0.1116 and
0.8042, respectively. The median hazard ratio is 1.0675. For an adverse
event, a hazard ratio below 1 favours treatment; the corresponding
no-effect value on the log scale is zero. The two rows transform both
the estimates and credible limits consistently.

``` r

ggplot(data.frame(hazard_ratio = exp(log_hr)), aes(hazard_ratio)) +
  geom_density(colour = "#0072B2", linewidth = 0.9) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#666666") +
  theme_bw() + labs(x = "Hazard ratio (treatment / control)", y = "Density")
```

![Posterior hazard ratio from the approximate common-effect calculation.
The dashed line marks a hazard ratio of
1.](bdpsurvival-vignette_files/figure-html/hazard-ratio-1.png)

Posterior hazard ratio from the approximate common-effect calculation.
The dashed line marks a hazard ratio of 1.

Repeat the two-arm analysis using draw-specific borrowing, retaining the
same settings. Treatment and control have separate weight distributions.

``` r

set.seed(42)
fit_two_mc <- bdpsurvival(Surv(time, status) ~ treatment,
  data = current_two, data0 = historical_two, breaks = c(3, 6, 12),
  a0 = 0.1, b0 = 0.1, method = "mc", fix_alpha = FALSE,
  discount_function = "identity", alpha_max = 1)
methods_two <- list("Fixed" = fit_two, "Draw-specific" = fit_two_mc)
knitr::kable(do.call(rbind, lapply(names(methods_two), function(method) {
  fit <- methods_two[[method]]
  hr <- exp(fit$final$posterior_loghazard)
  data.frame(Method = method,
    `Median treatment weight` = median(fit$posterior_treatment$alpha_discount),
    `Median control weight` = median(fit$posterior_control$alpha_discount),
    `Median HR` = median(hr), `95% CrI for HR` = format_interval(hr), check.names = FALSE)
})), digits = 3)
```

| Method | Median treatment weight | Median control weight | Median HR | 95% CrI for HR |
|:---|---:|---:|---:|:---|
| Fixed | 0.112 | 0.804 | 1.068 | 0.721 to 1.548 |
| Draw-specific | 0.099 | 0.487 | 1.006 | 0.667 to 1.563 |

## Piecewise exponential model

Partition follow-up at 0=\tau_0\<\tau_1\<\cdots\<\tau_J=\infty. For the
event-counting convention used when splitting records, interval j is
(\tau\_{j-1},\tau_j\]. Let d\_{ij} indicate an event in that interval
and t\_{ij} be exposure *within* the interval, not total follow-up.
Define D_j=\sum_i d\_{ij} and T_j=\sum_i t\_{ij}.

The hazard is constant at \lambda_j in interval j. Independent Gamma
priors with shape a_0 and rate b_0 give

\lambda_j\mid D_j,T_j\sim\operatorname{Gamma}(a_0+D_j,b_0+T_j).

The defaults are a_0=b_0=0.1. This prior is not flat. The rate parameter
has exposure-time units, so a change from months to days requires a
corresponding change in b_0 to preserve the same prior information.

The prior shape a_0 (Latin “a”) and the historical weight \alpha (Greek
“alpha”) have different roles:

| Quantity | R name | Role |
|:---|:---|:---|
| a_0 | `a0` | Baseline Gamma prior shape; default 0.1 |
| b_0 | `b0` | Baseline Gamma prior rate; default 0.1 |
| \alpha | `alpha_discount` in the fitted arm | Historical borrowing weight |
| \alpha\_{\max} | `alpha_max` | Maximum permitted weight, or the prescribed weight when `fix_alpha = TRUE`; default 1 |

Let \ell_j(t)=\max\\0,\min(t,\tau_j)-\tau\_{j-1}\\ be the length of
interval j before t. Then

S(t)=\exp\left\\-\sum_j\lambda_j\ell_j(t)\right\\.

The exported helper `ppexp(t, hazards, cuts = c(0, breaks))` returns
1-S(t) for each row of hazard draws. The fitting function calculates
survival directly from cumulative hazards. Borrowing comparisons use
retained log hazards, so they can distinguish draws even when the
corresponding survival probabilities round to zero or one.

### Comparing current and historical data

Comparisons use independent draws from the separate current and
historical posteriors within an arm. Write A_j=a_0+D_j and B_j=b_0+T_j,
with A\_{0j},B\_{0j} defined analogously for historical data. Gamma
shape-rate moments give \operatorname{Var}(\lambda_j)=A_j/B_j^2 and
\operatorname{Var}(\log\lambda_j)=\psi_1(A_j), where \psi_1 is the
trigamma function (Wellner 2008). The latter is exact; replacing it by
1/A_j is unsuitable for small shapes.

#### One arm: survival at the chosen horizon

At horizon t^\*, the log survival difference and its posterior variance
are

\Delta^{(s)}=\log S^{(s)}(t^\*)-\log S_0^{(s)}(t^\*)
=\sum_j\ell_j(t^\*)\\\lambda\_{0j}^{(s)}-\lambda_j^{(s)}\\,

V=\sum_j\ell_j(t^\*)^2\left(\frac{A_j}{B_j^2}+\frac{A\_{0j}}{B\_{0j}^2}\right),
\qquad Z^{(s)}=\frac{\Delta^{(s)}}{\sqrt{V}}.

Only intervals before the horizon contribute. For example, a horizon at
a cutpoint excludes the subsequent interval. Changing the horizon can
change the comparison, although the resulting historical weight is
applied to all intervals.

#### Two arms: compare studies separately within each arm

Within treatment, and separately within control, define the interval log
hazard differences and variances

C_j^{(s)}=\log\lambda\_{0j}^{(s)}-\log\lambda_j^{(s)}, \qquad
v_j=\psi_1(A\_{0j})+\psi_1(A_j).

Combine intervals using inverse-variance weights:

C^{(s)}=\frac{\sum_j C_j^{(s)}/v_j}{\sum_j1/v_j}, \qquad
V=\frac{1}{\sum_j1/v_j}, \qquad Z^{(s)}=\frac{C^{(s)}}{\sqrt{V}}.

Both borrowing methods use these analytic interval weights. Opposing
current–historical differences across intervals can cancel in the
combined comparison; inspect interval-specific hazards rather than
treating a high weight as evidence of agreement throughout follow-up.

#### Transforming the comparison into a weight

For `method = "fixed"`, estimate the directional probability
q=\Pr\\\Delta\>0\mid\text{data}\\ for one arm, or
q=\Pr\\C\>0\mid\text{data}\\ within each arm of a two-arm analysis. Then
use

p=2\min(q,1-q),\qquad \alpha=\alpha\_{\max}W(p).

For `method = "mc"`, use each standardized draw directly:

p^{(s)}=2\\1-\Phi(\|Z^{(s)}\|)\\, \qquad
\alpha^{(s)}=\alpha\_{\max}W(p^{(s)}).

The normal-tail transformation is applied once: Z^{(s)}=0 gives
p^{(s)}=1, and larger absolute differences give smaller values. It
remains an approximate comparison score, particularly with sparse,
skewed Gamma posteriors; exact variance standardization does not make
those distributions normal. Identical current and historical data
produce independent posterior draws, so their draw-specific weights need
not all equal one.

Neither transformation is a posterior probability of equality or
equivalence. `p_hat` contains the selected method’s transformed
comparison, including when `fix_alpha = TRUE`. In that case, the weight
is prescribed as `alpha_max` and does not depend on `p_hat`.

Log-Gamma sampling avoids an artificial positive hazard offset. Very
small hazards may still round to zero in the returned ordinary-scale
matrices, while their retained log values remain available for the
comparisons and treatment effect. More simulation draws reduce Monte
Carlo error; they do not remove the spread of draw-specific weights or
resolve the normal approximation’s limitations.

### Historical likelihood augmentation

Let L_j(\lambda_j) and L\_{0j}(\lambda_j) be the current and historical
interval likelihoods. For a chosen weight \alpha, historical likelihood
discounting combines them with the baseline prior \pi_0 as (Ibrahim and
Chen 2000)

\pi(\lambda_j\mid\text{data},\alpha) \propto
\pi_0(\lambda_j)L_j(\lambda_j)L\_{0j}(\lambda_j)^\alpha.

The historical likelihood contributes \lambda_j^{\alpha
D\_{0j}}\exp(-\alpha T\_{0j}\lambda_j). With the Gamma shape-rate prior,
the resulting update is

\lambda_j\mid\text{data},\alpha\sim\operatorname{Gamma}\left(
a_0+D_j+\alpha D\_{0j}, b_0+T_j+\alpha T\_{0j}\right).

The baseline prior enters once. The two separate source posteriors used
to estimate borrowing weights do not contribute additional copies of
that prior to augmentation. At \alpha=0, the update is the current-only
posterior. At \alpha=1, it is the pooled posterior under a common hazard
in each interval, with the same baseline prior. In two arms, treatment
and control use their own weights and hazard parameters.

With `method = "mc"`, substitute \alpha^{(s)} for each draw. The
resulting marginal distribution mixes over the generated weights and is
generally not one Gamma distribution. Its variance includes both the
conditional Gamma variance and variation between conditional means;
substituting the mean weight into one Gamma distribution does not
reproduce that mixture.

### Approximate two-arm effect

From the augmented arm hazards, calculate for each draw

R_j^{(s)}=\log\lambda\_{jT}^{(s)}-\log\lambda\_{jC}^{(s)}, \qquad
V_j=\operatorname{Var}\_s(R_j^{(s)}),

and combine intervals as

\beta^{(s)}=\frac{\sum_j R_j^{(s)}/V_j}{\sum_j1/V_j}.

These empirical variances summarize the *augmented* draws, including
variation from draw-specific weights. They serve a different purpose
from the analytic source-posterior variances used to choose borrowing
weights above.

This is the implemented inverse-variance approximation to a common log
hazard ratio. It is not an exact joint posterior fit of a
proportional-hazards Poisson regression. A common treatment effect
across intervals is needed for a common-hazard-ratio interpretation;
with time-varying effects it is a precision-weighted summary.
Large-event approximations such as 1/D\_{jT}+1/D\_{jC} do not describe
the implemented empirical variances exactly and are unsuitable when an
interval has no events.

## Assumptions and sensitivity

Assume independent participants and censoring independent of the event
process within each modeled group. Use the same event definition, time
origin, and time units across studies. Inspect follow-up support at the
chosen horizon; extrapolation beyond observed follow-up depends on the
piecewise hazard assumptions.

Check sensitivity to interval boundaries, survival horizon, Gamma
priors, borrowing method, discount function, and maximum weight. Sparse
or zero-event intervals can be strongly prior-dependent. Compare
interval-specific effects before interpreting a common hazard ratio, and
check numerical precision with more draws and independent seeds.

Trial-design use requires evaluation of the complete procedure,
including its operating characteristics under agreement and conflict
between studies. These worked examples do not provide that validation.

## References

Ibrahim, Joseph G., and Ming-Hui Chen. 2000. “Power Prior Distributions
for Regression Models.” *Statistical Science* 15 (1): 46–60.
<https://doi.org/10.1214/ss/1009212673>.

Wellner, Jon A. 2008. *Gamma, Digamma, and Polygamma*. University of
Washington, Statistics 581 handout.
<https://sites.stat.washington.edu/jaw/COURSES/580s/581/HO/digamma2.pdf>.
