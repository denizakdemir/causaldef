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
  implementing stabilized inverse probability weighted Cox models.

### New Identification Methods

* **Front-Door Kernel** (`frontdoor_effect()`): Implements Theorem 2.2 with 
  plugin and doubly-robust estimators for front-door identification.
  
* **Transport Deficiency** (`transport_deficiency()`): Measures distribution 
  shift between source and target populations with IPTW and calibration methods.

* **Instrumental Variables** (`iv_effect()`): Full IV support with 2SLS, Wald,
  and LIML estimators. Includes weak instrument diagnostics and validity tests
  via `test_instrument()`.

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
  diagnostics with Theorem 5.2 implementation.
  
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
