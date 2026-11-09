## Submission type
This is an update from the current CRAN version, 0.2.1, to 0.3.0.

## Test environments
* local Ubuntu 20.04 (x86_64-pc-linux-gnu), R 3.6.3
* (to add before submission: win-builder R-devel and a current R 4.x run)

## R CMD check results
`R CMD check --as-cran` with `_R_CHECK_FORCE_SUGGESTS_=false`

0 errors | 1 warning | 2 notes, all from the local environment:

* WARNING: `qpdf` is not installed locally (needed only for the PDF size-reduction check).
* NOTE: Suggested packages not available locally: `MatchIt`, `tmle`, `SuperLearner`, `grf`, `plumber`, `cmprsk`.
* NOTE: "unable to verify current time" (no network time source in the local check).

The slow `confounding_frontier()` example (previously over 5 s) has been reduced to about 2.5 s.

## What changed since 0.2.1
This release corrects decision-theoretic constants that a re-check of the companion
manuscript showed to be wrong. Numerical outputs of the affected functions change:

* `policy_regret_bound()`: the transfer penalty is 2·|A|·M·δ (previously M·δ) and the minimax floor is (M/4)·δ (previously (M/2)·δ). New arguments `n_actions` and `policy_type`.
* `policy_regret_bound_vc()`: corrected finite-sample bound; the new argument `eta` is required.
* `rkhs_rate_bound()`: corrected rate (n·η)^(−β/(2β+d_W+d_Y)); new argument `d_y`.
* The manuscript is cited by its concept DOI (10.5281/zenodo.18367346).

See `NEWS.md` for details.

## Downstream dependencies
There are no downstream dependencies.
