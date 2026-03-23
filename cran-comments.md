## Test environments
* local linux (Ubuntu 20.04.6 LTS), R 3.6.3
* win-builder (Windows Server 2022 x64, R 4.5.3 ucrt)

## R CMD check results
`R CMD check --as-cran --no-manual` with `_R_CHECK_FORCE_SUGGESTS_=false`

0 errors | 1 warning | 2 notes

win-builder (R-release):

0 errors | 0 warnings | 1 note

Notes:
* This is a new submission.
* Suggested packages not available in the local check environment:
  `MatchIt`, `tmle`, `SuperLearner`, `grf`, `plumber`, `cmprsk`
* win-builder also reports DESCRIPTION spellcheck suggestions for
  `Akdemir` and `computable`, which are not package issues.

Warning:
* `qpdf` is not installed in the local environment, so PDF size reduction checks
  could not be run.

## Responses to CRAN Feedback
* Clarified the package contract to distinguish theorem-backed bounds from
  computable deficiency proxies, sensitivity diagnostics, and heuristic modules.
* Disabled previously misleading estimator labels:
  `frontdoor_effect(method = "dr")` and `iv_effect(method = "liml")`.
* Reworked `nc_diagnostic()` to use permutation-based screening and to report
  the `kappa` sensitivity bound separately.
* Added explicit survival runtime guards on older R versions where the current
  `survival` package fails internally because `base::deparse1` is unavailable.
* Updated vignettes and examples so survival code is evaluated only on
  compatible runtimes.

## Downstream dependencies
There are no downstream dependencies.
