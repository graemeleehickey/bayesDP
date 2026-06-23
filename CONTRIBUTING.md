# Contributing to bayesDP

Thanks for your interest in contributing to bayesDP! This document
describes a few conventions to follow when submitting changes.

## Pull requests

- Please raise an issue describing the bug or feature before opening a
  large pull request, so the approach can be discussed first.
- Keep pull requests focused: one logical change per pull request makes
  review easier.
- Run `R CMD check` (e.g. via `devtools::check()`) locally and ensure
  there are no new errors, warnings, or notes before submitting.
- Add or update tests under `tests/testthat/` for any behaviour you
  change. Where a result is analytically tractable, pin computed values
  against the closed-form answer rather than only checking that the code
  runs.

## Code style

- Code is formatted following the [tidyverse style
  guide](https://style.tidyverse.org/). Running `styler::style_pkg()`
  before committing is encouraged.
- Documentation is written with roxygen2; edit the roxygen comments in
  the `R/` source files (not the generated `man/*.Rd` or `NAMESPACE`)
  and run `devtools::document()` to regenerate.

## NEWS.md

User-facing changes must be recorded in `NEWS.md` under the current
development version heading (the top-most `# bayesDP x.y.z.9000`
section).

Within each version, entries are grouped under level-two (`##`)
subsection headings so the changelog is easy to scan. Use the following
sections, in this order, and include only those that apply:

- **`## Bug fixes`** — corrections to incorrect behaviour.
- **`## New features`** — new functions, arguments, or capabilities.
- **`## Tests`** — additions or changes to the test suite and coverage.
- **`## Housekeeping`** — refactoring, dependency and CI updates,
  documentation tidy-ups, deprecation fixes, and other internal changes
  with no direct effect on results.

Each entry is a single bullet point. Reference the relevant GitHub issue
or pull request number where helpful. For example:

``` markdown
# bayesDP 1.3.7.9000

## Bug fixes

* Fixed `bdplogit()` failing to locate the response in the model fit (#12).

## Housekeeping

* Updated GitHub Actions workflows to the latest `r-lib/actions` examples.
```

Older released versions in `NEWS.md` predate this convention and are
left as-is; the subsection structure applies going forward.
