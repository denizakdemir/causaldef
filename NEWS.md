# causaldef 0.3.0

## Corrections to decision-theoretic bounds (numerical outputs change)

A full re-check of the companion manuscript found that two constants used by the
policy-bound utilities were wrong; counterexamples are in the manuscript (Zenodo
concept DOI 10.5281/zenodo.18367346, version 5).

* **`policy_regret_bound()`**: the transfer penalty is now
  2·|A|·M·δ (4Mδ for a binary treatment with covariate-dependent
  policies; 2Mδ when comparing fixed actions via `policy_type = "constant"`). Versions
  up to 0.2.1 reported Mδ, which does not bound the regret: with two actions and no
  covariates the worst-case regret can be 2Mδ while the simulated regret is 0. The
  minimax floor is now (M/4)δ (previously (M/2)δ, which is false; 1/4 is sharp). New
  arguments `n_actions` and `policy_type`; new fields `value_error`,
  `transfer_multiplier` and `floor_multiplier`. `obs_regret` is the regret computed
  under the adjusted (simulated) interventional laws.
* **`policy_regret_bound_vc()`**: now implements the corrected finite-sample bound for
  inverse-propensity-weighted empirical welfare maximization with a known propensity
  in [η, 1 − η]: 4Mδ + C (M/η)(√(v/n) + √(log(1/ξ)/n)). The new argument `eta` is
  required. `obs_regret` is ignored, with a warning.
* **`rkhs_rate_bound()`**: now returns the rate (nη)^(−β/(2β + d_W + d_Y)), with the new
  argument `d_y` (0 for finite outcome spaces). The previous formula ignored the
  outcome dimension and misplaced η. `xi` is retained but unused.
* **`plot.policy_bound()`** (`type = "safety_curve"`) and the plumber `/policy-bound`
  endpoint use the corrected constants.
* **`confounding_frontier()` documentation**: the closed-form lower bound now shows the
  factor 1/2 from the two-point method. The computed frontier was already correct.

## Other

* The `confounding_frontier()` examples use smaller grids, so they run in about 2.5 s instead of 16 s.
* The manuscript is cited by its concept DOI 10.5281/zenodo.18367346, which always
  resolves to the latest version.

# causaldef 0.2.1

## Bug Fixes

* **`confounding_frontier()` / `.tv_distance_normal()`**: Fixed a numerical bug
  in the closed-form crossing-point TV computation that returned values near 0
  instead of near 1 when both Gaussian components have very small standard
  deviations relative to their mean separation (a near-degenerate point-mass
  regime where the quadratic discriminant loses precision). Replaced with
  integration split into pieces sized to each component's own sd. Verified
  against an independently-constructed cross-check across 3000 random cases
  spanning 8 orders of magnitude in sd (max discrepancy ~4e-8); the function's
  default usage range within `confounding_frontier()` was never affected by
  this bug.

## Documentation

* Removed an unspecified "universal constant" phrase from `rkhs_rate_bound()`
  and `confounding_frontier()` documentation (the underlying manuscript bound
  was found to not hold with any single universal constant for large `sigma_A`
  and was replaced with an explicit, unconditional bound).
* Hedged interpretive guidance in the README and vignettes.
* Reworded documentation, README, and vignettes to reduce narrative
  over-attribution to "Le Cam" (e.g. "Le Cam deficiency" to "deficiency",
  "Le Cam's theory" to "the classical theory of statistical experiment
  comparison"), while keeping the Le Cam & Yang (2000) bibliography citation.
  No functional code changes from this item.

# causaldef 0.2.0

## Major Changes

### New Estimation Methods

* **TMLE (Targeted Maximum Likelihood Estimation)**: Added `method = "tmle"` 
  for doubly-robust estimation using the `tmle` package with SuperLearner.
  
* **MatchIt Integration**: Proper propensity score matching via `method = "matching"` 
  using the `MatchIt` package with nearest neighbor matching.
  
* **Causal Forests (grf)**: Added `method = "grf"` for heterogeneous treatment 
  effect estimation using Generalized Random Forests.
  
* **Cox IPTW for Survival**: Added `method = "cox_iptw"` for survival outcomes, 
  implementing stabilized inverse probability weighted Cox models on
  compatible survival runtimes.

### New Identification Methods

* **Front-Door Kernel** (`frontdoor_effect()`): Implements the front-door kernel existence result (`thm:frontdoor`) with 
  a plugin estimator and a heuristic front-door deficiency proxy.
  
* **Transport Deficiency** (`transport_deficiency()`): Measures distribution 
  shift between source and target populations with proxy diagnostics.

* **Instrumental Variables** (`iv_effect()`): IV support with 2SLS and Wald
  estimators, plus weak instrument diagnostics and validity tests via
  `test_instrument()`.

### New Outcome Types

* **Competing Risks** (`causal_spec_competing()`): Full support for time-to-event
  data with multiple event types. Implements cause-specific and subdistribution
  hazard estimation via `estimate_deficiency_competing()`.

### Performance Improvements

* **Parallel Bootstrap**: New `parallel = TRUE` argument in `estimate_deficiency()` 
  enables parallel processing via `future.apply` for faster inference with 
  large bootstrap samples.

* **Stabilized IPTW Weights**: Propensity scores are now bounded to [0.01, 0.99] 
  to prevent extreme weights.

### Commercial Features

* **Shiny Dashboard** (`run_causaldef_app()`): Interactive web application for
  deficiency analysis with data upload, method comparison, and report export.
  
* **Standalone Deployment** (`create_shiny_app_files()`): Generate app files
  for shinyapps.io or Shiny Server deployment.

* **REST API** (`create_plumber_api()`, `run_causaldef_api()`): Full REST API
  via plumber for SaaS deployment. Includes endpoints for deficiency estimation,
  policy bounds, confounding frontiers, and transport analysis. Docker-ready.

### New Vignettes

* `negative_controls.Rmd`: Comprehensive guide to using negative control 
  diagnostics with the negative control bound (`thm:nc_bound`).
  
* `policy_learning.Rmd`: Guide to safe policy learning with decision-theoretic 
  bounds and the safety floor concept.

### Infrastructure

* Added comprehensive demo script: `inst/examples/complete_demo.R`
* Added pkgdown configuration for documentation website
* Added new test suites for all new functions
* Expanded DESCRIPTION suggests to include: `tmle`, `MatchIt`, `grf`, 
  `SuperLearner`, `future.apply`, `shiny`, `cmprsk`, `plumber`, `jsonlite`

## Bug Fixes

* Fixed weight normalization in `matching` method
* Improved fallback handling when optional packages fail
* Clarified package semantics so theorem-backed bounds, computable proxies,
  sensitivity diagnostics, and heuristic modules are labeled distinctly
* Disabled unsupported `frontdoor(method = "dr")` and `iv_effect(method = "liml")`
* Reworked `nc_diagnostic()` into permutation-based screening plus
  `kappa`-sensitivity bounds
* Added explicit survival runtime guards for older R versions where
  `survival` internals require `base::deparse1`
* Made `policy_regret_bound()` method selection explicit and recorded
  optimistic post-selection when used

---

# causaldef 0.1.3

* Initial CRAN submission version
* Core functions: `causal_spec()`, `causal_spec_survival()`, 
  `estimate_deficiency()`, `nc_diagnostic()`, `confounding_frontier()`, 
  `policy_regret_bound()`
* Basic methods: `unadjusted`, `iptw`, `aipw`
* Four vignettes covering methodology and survival analysis

---

# causaldef 0.1.0

* Initial development version
* Implements Le Cam deficiency theory for causal inference
* Based on Akdemir (2026) "Constraints on Causal Inference as Experiment Comparison"
