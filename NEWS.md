# bayesDP 1.3.9

## Bug fixes

* Fixed direct use of fitted models with `modelsummary()`: tidiers now accept its neutral `vcov = NULL` and `coef_rename = FALSE` arguments while retaining posterior SDs and equal-tailed credible intervals. Unsupported covariance overrides still produce an error; table labels can be changed with `coef_map`.
* Fixed treatment-only regression formulas: `bdplm()` no longer adds unintended covariate columns under `method = "mc"`, and logistic posterior draws contain only the intercept and treatment effect. Posterior columns are selected by name in both regressions.
* `bdplogit()` now defaults to its supported `"fixed"` method and uses distinct seeds drawn from R's random-number state for its two comparisons and final fit. `set.seed()` controls reproducibility, and stage seeds are retained in `args1$mcmc_seeds`. Results will differ from earlier versions that always used MCMCpack's default seed.
* `bdpsurvival()` now uses evaluated survival responses, supporting transformed times, event expressions, named arguments, and `survival::Surv()`. Missing responses, unsupported censoring types, and invalid treatment codes cause explicit errors before interval splitting. Valid one-arm and historical-only-control analyses remain supported.
* Added validation of binomial and normal summaries and borrowing settings before sampling. Impossible counts, invalid sample sizes or SDs, non-finite values, incomplete groups, and unintended vector inputs now cause errors.
* `bdplogit()` now stops with a stage-specific explanation when complete or quasi-complete separation, or a rank-deficient model matrix, prevents a supported analysis. Checks cover each arm's borrowing comparison, the historical regression when arm priors are estimated, and the current-data regression. This prevents borrowing weights from improper comparison posteriors and historical priors derived from infinite maximum likelihood estimates. The current-data check reflects a limitation of the sampler, including when proper priors are supplied.
* With `fix_alpha = TRUE`, `bdplogit()` now uses the prescribed arm weights without fitting comparison models; internal comparison probabilities are recorded as missing. Historical-prior estimation and current-data checks still apply. Prescribed weights must be finite values in [0, 1].
* `bdplm()` and `bdplogit()` now reject formula offsets before fitting instead of silently omitting them. This includes `stats::offset()`, zero offsets, and cancelling offsets, regardless of the borrowing settings or supplied historical priors. Offset models fitted with earlier versions must be reconsidered; their results do not represent the specified model.
* Corrected survival augmentation to discount historical event counts and exposure while including the baseline Gamma prior once (#17). The augmented shape is now `a0 + D + alpha * D0`, and the rate is `b0 + T + alpha * T0`, removing the extra `alpha * a0` and `alpha * b0` contributions.
* Applied the correction to one- and two-arm analyses, including prescribed, estimated fixed, and draw-specific weights. Full borrowing now gives the pooled interval-hazard posterior with the same baseline prior and interval boundaries. Survival estimates and treatment effects may change, especially with sparse events or influential priors; refit earlier analyses to use the corrected update.
* Corrected survival Monte Carlo borrowing to standardize log survival differences by their posterior SDs and apply the two-sided normal-tail transformation once (#16).
* Current–historical comparisons within each arm now use analytic log-Gamma variances and inverse-variance interval weights. This changes two-arm borrowing under both `method = "fixed"` and `method = "mc"` (#16).
* Removed the artificial hazard offset and retained log-Gamma draws for stable comparisons and two-arm effects with sparse or zero-event intervals. One-arm comparisons use cumulative hazards directly, avoiding loss of information when survival probabilities round to zero or one. The changed sampling calculations also change results obtained with a previous random seed.
* Made `p_hat` consistent with the chosen comparison method when `fix_alpha = TRUE`, while preserving the prescribed weight. Added checks for invalid prior, horizon, and borrowing settings, and corrected interval handling for a single constant hazard and intervals with no observations.
* Corrected survival summaries to count participants and calculate risk sets from their observed follow-up, rather than treating split hazard intervals as separate participants (#18). Artificial interval endpoints no longer appear as observation times in the one-arm summary.
* Corrected regression residual summaries to use covariates and the reported intercept on the same original scale, including encoded factor covariates (#18).
* Applied the correction to logistic regression, which shares the regression summary method. Its residuals now use fitted probabilities, all coefficient rows are retained, and no residual SD is reported for binary outcomes.
* Survival fits now retain participant-level records, and regression fits retain covariate centering means. Older saved fits must be refitted before using the affected summaries; the functions explain when this is necessary.
* Corrected the normal-data `probability_discount()` helper to standardize Monte Carlo mean differences by their conditional standard error, matching `bdpnormal()` and preserving the comparison when outcome units change.
* Corrected discount plots to include `alpha_max` for each arm and to show a constant weight when `fix_alpha = TRUE`.
* Corrected the fixed-weight probability definitions and the normal vignette's treatment-effect conclusion (#19, #20).
* Fixed equation rendering in the pkgdown articles with KaTeX and updated the standalone HTML vignettes to MathJax 3.
* Corrected the diagonal-matrix notation in the linear regression vignette so that it renders correctly on the website.

## New features

* Added `centering = "pooled"`, `"current"`, or `"historical"` to `bdplm()` and `bdplogit()` to choose which observations define the covariate means. The default preserves pooled centering. Each choice uses both treatment arms, without discount weighting, and applies the same means to both datasets. Fits retain the choice and means; changing the reference can change the posterior under the independent historical arm priors.

* Added `coef()`, `vcov()`, and `confint()` for linear and logistic regression fits. These return posterior mean coefficients, their posterior covariance, and equal-tailed credible intervals on the original covariate scale. Added `nobs()` for all fitted models, counting current participants separately from historical participants.
* Added `posterior::as_draws_df()` and `posterior::as_draws_matrix()` for all fitted models. Conversions preserve retained draws and their order, including paired arm differences, regression residual SDs, and survival interval log hazards. The posterior package is optional.
* Added `tidy()` and `glance()` methods for all fitted models, also available through broom. `tidy()` reports posterior estimates, SDs, and optional equal-tailed credible limits; logistic and two-arm survival effects can be exponentiated. `glance()` reports current and historical sample sizes, retained draw counts, and arm-specific discount weights, with event counts for survival fits.
* Numerical summaries are shared with `summary()`. The new methods use saved draws without refitting or changing the random-number state. Participant-level `augment()` methods are not provided.

* Added `type = "effect"` plots for binomial and normal treatment-minus-control contrasts. Posterior plots now use explicit parameter labels and consistent information-source colours and line types.
* Added a light/dark mode switch to the pkgdown website, including matching colours for R code.

## Tests

* Added checks for coefficient extraction, joint posterior covariance, credible intervals, participant counts, posterior draw conversion, and direct modelsummary tables across the supported model types.
* Added tests for treatment-only models, logistic seed control and method defaults, equivalent survival-response expressions, preservation of participant counts, invalid treatment coding, and summary-data input boundaries.
* Expanded checks for incomplete inputs, historical-only controls, factor contrasts, regression discount functions, plot options, and survival calculations at interval boundaries and numerical limits. Added comparisons with analytical posterior moments and integrated hazards.
* Added logistic regression checks for zero-event and all-event groups, covariate separation, quasi-complete separation, rank deficiency, prescribed weights, and supplied historical priors.
* Added checks that regression offsets are rejected before sampling across supported fitting settings and that ordinary covariates named `offset` remain supported.
* Added checks against analytical Gamma and survival moments, current-only and pooled analyses, different arm weights, draw-specific Gamma mixtures, zero events, and intervals without observations. The tests use influential priors as well as the default prior, so duplicated prior contributions are detected.
* Added independent calculation checks, distributional checks for log-Gamma draws, and tests covering time-unit changes, interval and draw ordering, zero events, extreme survival probabilities, all discount functions, and prescribed weights.
* Added tests for participant counts, risk sets, tied and censored observations, events at interval boundaries, historical-only controls, and regression residuals with shifted and factor covariates.
* Added checks against numerical integration, outcome-unit changes, the model's comparison draws, discount-plot limits, treatment-effect plots, numerical summary labels, and the documented regression covariate-prior behaviour.
* Included the test suite in source packages so `R CMD check` runs it, and declared the existing `withr` test dependency.

## Housekeeping

* Added `detectseparation` for logistic separation checks and documented the supported fitting conditions and the distinction between prescribed and estimated weights. Logistic sampling still requires `method = "fixed"`.
* Updated the survival and borrowing vignettes to explain historical likelihood discounting, distinguish the historical weight `alpha` from the Gamma shape `a0`, and interpret current-only and full-borrowing limits. Regenerated the worked examples with the corrected update.
* Updated the survival vignette with the corrected formulas, a comparison of fixed and draw-specific borrowing, numerical considerations, and the limitations of normal-tail comparisons with sparse events.
* Updated the regression and survival vignettes to demonstrate the corrected summaries, and labelled one-arm survival uncertainty as posterior SDs and equal-tailed credible intervals.
* Reorganised the four model vignettes around worked examples, interpretation, assumptions, and model details, with a shared historical-borrowing overview and formal references.
* Added comparisons of current-only, estimated, and full borrowing; generated numerical interpretations from fitted draws; and included a survival example with right censoring.
* Clarified the distinction between estimated fixed and draw-specific weights, conditional and marginal distributions, baseline priors, regression covariate-prior precision, and posterior medians versus means.
* Labelled binomial and normal summary intervals as equal-tailed credible intervals and retained four decimal places for arm estimates.
* Added strict KaTeX validation to the website build so unsupported LaTeX commands and invalid equations prevent deployment.
* Improved equation spacing and alignment, split long equations across lines, and divided dense paragraphs in all four vignettes to make them easier to read.
* Standardised spacing in R examples and documented function calls throughout the vignettes, README, and help pages.
* Added the CRAN DOI badge to the README and package website.

# bayesDP 1.3.8

## Bug fixes

* Fixed `bdplm()` and `bdplogit()` producing invalid historical borrowing when covariates were not mean-centered. Because the models use an intercept-free parameterization (separate treatment and control means), uncentered covariates made the arm-mean estimators strongly correlated and inflated their standard errors as extrapolation errors at covariate = 0, corrupting the (diagonal) discount prior. Both functions now automatically mean-center covariates on their pooled (current plus historical) mean and back-transform the reported intercept, so estimates are invariant to covariate location shifts (#1).
* Fixed a bug where the `summary` methods for `bdpnormal` and `bdpbinomial` read the one-/two-arm indicator from the wrong list element (`args` instead of `args1`), causing two-arm fits to be summarised as one-arm.
* Fixed a bug in the two-arm `bdpsurvival` summary that errored when current control data were absent.
* Fixed an invalid matrix index (`Y[, 0]`) used when computing the default `surv_time` in `bdpsurvival`.
* Fixed the `compare` argument being silently dropped (passed into `paste()`) rather than stored in the `bdpnormal` and `bdpbinomial` fit objects.
* Fixed the `mc` discount-weight Z-statistic in `bdplm` to divide by the standard error rather than the variance.
* Fixed `plot` methods hanging on the interactive "Hit <Return>" prompt in non-interactive sessions (e.g. tests, CI); `par(ask = ...)` now respects `interactive()`.
* Fixed `bdplogit()` failing during its main model fit because the analysis data passed to `MCMClogit()` omitted the response variable. The discounted prior precision matrix is now also passed to `MCMClogit()` correctly (#12).
* Fixed `alpha_discount()` so `alpha_max` is respected when `discount_function = "identity"` (#6).
* Fixed `bdpnormal()` one-arm normal fits with only one source of data for an arm (current-only or historical-only internally) returning an over-dispersed `posterior_mu`. These branches now return the conjugate posterior of the mean rather than adding an extra observation-level draw (posterior-predictive-like variance) (#15).

## Tests

* Re-enabled the `testthat` test harness (`tests/testthat.R`).
* Rewrote the test suite with `expect_*()` assertions: augmented one-arm binomial and normal posterior means are pinned against their closed-form conjugate values, and the fixed bugs (one-/two-arm dispatch, stored `compare` flag, default survival time, two-arm survival summary) are now guarded by tests. Plot calls in tests pass an explicit `type` so they no longer prompt for input.
* Expanded test coverage from ~60% to ~76%, adding tests for `alpha_discount()` and `probability_discount()` (both now fully covered), the `ppexp()` vector and matrix paths, the `print` methods (now fully covered), additional `plot` branches, input-validation paths, the `bdplogit()` main fit path, factor-covariate handling in `bdplm()` and `bdplogit()`, and the `mc` discounting method for `bdpnormal` and `bdpbinomial`.
* Added regression tests pinning the `bdpnormal` flat-prior draw of the mean (`posterior_flat_mu`) and the fixed current-only `posterior_mu` against their closed-form conjugate (Student-t) variance.
* Expanded `method="mc"` documentation in the binomial, normal, and survival interfaces/vignettes to note that per-iteration recomputation of the stochastic comparison yields a random `alpha_discount` sequence that can show noticeable Monte Carlo variability (#4).
* Guarded the plotting tests with a null graphics device so they no longer write a stray `Rplots.pdf`.

## Documentation

* Clarified the `prior_covariate_sd` documentation in `bdplm()` and `bdplogit()` to note that covariate effects carry an intentional near-zero discount weight, making their priors effectively flat. The supplied value has negligible influence on the posterior, and the effective prior standard deviation is roughly 1e6 larger than the nominal value at the default (#2).

## Housekeeping

* Removed AppVeyor continuous integration.
* Updated GitHub Actions workflows to the latest `r-lib/actions` examples (`actions/checkout@v6`, `codecov/codecov-action@v6`, `actions/upload-artifact@v7`).
* Made Codecov upload issues non-fatal in the coverage workflow so tests and coverage generation remain the CI gate while Codecov service/signature failures do not fail the build.
* Added a `pkgdown` website and accompanying GitHub Actions workflow.
* Replaced deprecated `ggplot2::aes_string()` with `aes()` and the `.data` pronoun in all `plot` methods.
* Replaced deprecated `size` aesthetic with `linewidth` in `geom_line()` calls.
* De-duplicated the internal `model.matrixBayes()` helper (previously defined identically in both `bdplm` and `bdplogit`) into a single internal file.
* Removed leftover commented-out debugging code.
* Collapsed a redundant conditional in `posterior_survival()` where both branches initialized identical hazard matrices (#8).
* Tidied the `mc` sigma2 sampling in `bdplm()` and kept the sampling weights aligned with the candidate grid when some marginal log-likelihoods are non-finite (#9).
* `ppexp()` now validates its `x` argument and errors with an informative message when it is neither a numeric vector nor a matrix (#10).
* Survival curves in the `plot` and `summary` methods are now computed with a vectorised C++ routine (`ppexpMV`) that transposes the hazard matrix once across all time points, instead of looping `ppexp()` per time point (#11).
* Avoided recomputing the per-interval sufficient statistics in `posterior_survival()`; the augmentation step now reuses the values already computed during the discount phase (#7).
* Removed a redundant `useDynLib()` directive in the package namespace.
* Added contributor guidance documenting the `NEWS.md` subsection convention for future releases.
* Expanded `README.md` with links, supported analyses, examples, and citation guidance (#13).
* Added the CRAN package URL to `DESCRIPTION` (#14).
* Clarified in `posterior_normal()` that the flat-prior draw of the mean is the conjugate posterior (scale `sqrt(sigma^2 / N)`), not the posterior predictive.
* Gave each vignette a descriptive title (previously all titled "BayesDP") and removed unused `params`/`EVAL` scaffolding from the vignette headers.

# bayesDP 1.3.7

* Updated GitHub actions workflows.
* Updated README badges.
* Fixed .Rd file itemize list for `bdpbinomial`, `bdpnormal`, and `bdpsurvival`.
* Fixed logical check in `bdpsurvival`.
* Fixed spelling mistakes in documentation.
* Minor formatting updates.
* Add reverse dependency checks.

# bayesDP 1.3.6

* Fixed CRAN CMD warnings for S4 generics.

# bayesDP 1.3.5

* Fixed CRAN CMD warnings for S4 generics.

# bayesDP 1.3.4

* Updates to README and DESCRIPTION.
* Updates to .gitignore file.
* Add codecov for coverage checking.
* Updates to R code formatting.
* Added GitHub Actions CI.

# bayesDP 1.3.3

* New package maintainer (Graeme L. Hickey) since package was orphaned.
* Updates to README, DESCRIPTION, NAMESPACE.
* Added `stop` break to `discount_logit` for `method = mc`.

# bayesDP 1.3.2

* Minor `bdplm` vignette typo fixes.

# bayesDP 1.3.1

## Major new features

* Changes to inputs for `bdpsurvival`.
  + Current and (optional) historical data are specified in separate data frames.
* Updated normal approximation used for `method = "mc"` of the `bdpbinomial` and `bdpnormal` functions.

## Bug fixes and minor improvements

* Summary method for `bdplm` now exists and mimics `lm`.
* Removed `bdpbinomial` vignette language around success (vs event).
* Reported one-arm sample size for `bdpsurvival` print method adjusted to current data only.

# bayesDP 1.3.0

## Major new features

* Addition of the `bdplm` function for two-arm trials.
* Users can now choose between 3 discount functions via the `discount_function` input.
  + Weibull CDF.
  + Scaled Weibull CDF - scales the Weibull CDF so that the max possible value is 1.
  + Identity - sets the discount weight to the posterior probability.
* Removal of `bdpregression`.

## Bug fixes and minor improvements

* Removed two-sided and one-sided function inputs to avoid confusion.
* Posterior probabilities for `method = "mc"` switched from `pshisq` to pnorm.
* Updated vignettes to reflect new features.

# bayesDP 1.2.0

## Major new features

* Supports one-arm regression analysis.
* Two additional modular functions.
* Implementation of Monte Carlo-based estimation of alpha discount.

## Bug fixes and minor improvements

* Fixes to class slots.
* Added `print` input to plot method.

# bayesDP 1.1.0

## Major new features

* Supports two-arm survival analysis via hazard rate comparisons.
* Completely revamped summary and print methods to produce better formatted results.
* Plot method allows users to specify a `type`.
* Added vignettes for each of `bdpbinomial`, `bdpnormal`, and `bdpsurvival`.
* Implemented the `fix_alpha` input which allows users to set the historical data weight at `alpha_max`.

## Bug fixes and minor improvements

* Fixed error with two-arm analysis where models did not fit if either the current or historical control data were not input.
* Changed `two_side` input to logical.
* Consolidated several internal functions into a single function for computational efficiency gains.

# bayesDP 1.0.3

* README update.
* Added plot types.
* Added Vignettes.
* Added logo.
* Improved documentation.
* Updated `print`, `summary`, `plot` methods.
* Refactored `bdpnormal` / `bdpbinomial`.

# bayesDP 1.0.2

* Crucial bugfixes.

# bayesDP 1.0.1

* User requested bugfixes.

# bayesDP 1.0.0

* Initial CRAN release with normal, binomial and survival functions.
