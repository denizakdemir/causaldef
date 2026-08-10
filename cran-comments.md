## Submission type
This is an update from the current CRAN version, 0.2.0, to 0.2.1.

## Test environments
* local macOS (aarch64-apple-darwin20), R 4.4.1

## R CMD check results
`R CMD check --as-cran` with `_R_CHECK_FORCE_SUGGESTS_=false`

0 errors | 1 warning | 2 notes

Notes:
* Suggested packages not available in the local check environment:
  `MatchIt`, `tmle`, `SuperLearner`, `grf`
* Local HTML manual check reports `<main>` as unrecognized by the installed
  `tidy`; this is a known local-toolchain limitation (outdated HTML5 support
  in system `tidy`), not a package issue.

Warning:
* CRAN incoming feasibility check locally compares against the last-known
  CRAN version and will clear once 0.2.1 is the version under review.

## What changed since 0.2.0
* Bug fix: `confounding_frontier()` / internal `.tv_distance_normal()`
  returned incorrect (near-0 instead of near-1) total-variation values in a
  near-degenerate regime where both Gaussian components have very small
  standard deviations relative to their mean separation. Fixed via a
  component-scaled split-integration approach; verified against an
  independent cross-check across 3000 random cases spanning 8 orders of
  magnitude in sd (max discrepancy ~4e-8). The function's default usage range
  within `confounding_frontier()` was never affected by this bug.
* Documentation: removed an unspecified "universal constant" phrase from
  `rkhs_rate_bound()` and `confounding_frontier()` docs (replaced with an
  explicit, unconditional bound); hedged interpretive guidance in the README
  and vignettes; reworded documentation to reduce narrative over-attribution
  to "Le Cam" while keeping the Le Cam & Yang (2000) citation. No functional
  code changes from the documentation items.

See `NEWS.md` for the full changelog.

## Downstream dependencies
There are no downstream dependencies.
