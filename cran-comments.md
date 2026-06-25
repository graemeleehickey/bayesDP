## Submission summary

This is a bug-fix release. It corrects several issues in the model-fitting
functions, most notably invalid historical borrowing in `bdplm()` and
`bdplogit()` when covariates were not mean-centered (these functions now
mean-center covariates internally and back-transform the reported intercept).
See NEWS.md for the full list of changes.

## Test environments

* local macOS (Tahoe 26.5), R 4.5.3
* ubuntu (via GitHub Actions, release + devel)
* macOS (via GitHub Actions, release)
* windows (via GitHub Actions, release)
* windows (via win-builder, old + release + devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies for this package on CRAN.
`revdepcheck::revdep_check()` checked 0 reverse dependencies.

 * We saw 0 new problems
 * We failed to check 0 packages
