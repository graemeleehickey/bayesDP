## Submission summary

This updates bayesDP from 1.3.8 to 1.3.9. It corrects survival historical borrowing and summaries, regression fitting and residual summaries, and input validation. It also adds a choice of covariate centering and methods for extracting coefficients, credible intervals, posterior draws, and tabular summaries. The vignettes and help files explain the corrected calculations and their interpretation.

Some corrections change numerical results, and affected analyses should be refitted as described in NEWS.md. The maintainer and licence are unchanged.

## Test environments

* Local macOS Tahoe 26.7.3 (arm64), R 4.5.3
* Ubuntu 26.04.1 LTS, R 4.6.2, R-devel (2026-11-23 r90971), and R 4.5.3, via GitHub Actions
* macOS Tahoe 26.7.2 (arm64), R 4.6.2, via GitHub Actions
* Windows Server 2022 (x64), R 4.6.2, via GitHub Actions

## R CMD check results

All five GitHub Actions checks completed with 0 errors, 0 warnings, and 0 notes, using `--as-cran --no-manual`.

The submission archive was built with R 4.6.2 and additionally checked locally using `--as-cran --run-donttest`, including PDF and HTML manual checks, CRAN incoming checks, and all suggested dependencies. This completed with 0 errors, 0 warnings, and 1 note:

> unable to verify current time

This note arose from R's external clock check. The timestamp check passed in all five GitHub Actions environments. All local package checks passed, including 2,639 test assertions with no failures, warnings, or skips. The documentation URL check also passed.

## Reverse dependencies

The CRAN package database was checked on 24 November 2026. No packages list bayesDP in Depends, Imports, LinkingTo, Suggests, or Enhances, so no reverse-dependency package checks were required.
