
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causaldef <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/denizakdemir/causaldef)
[![CRAN status](https://www.r-pkg.org/badges/version/causaldef)](https://CRAN.R-project.org/package=causaldef)
<!-- badges: end -->

**causaldef** implements Le Cam deficiency theory for causal inference, providing quantitative bounds on information loss from confounding, selection bias, and distributional shift.

Unlike traditional sensitivity analysis which focuses on *"how much bias"* exists, `causaldef` answers the decision-theoretic question: **"how much regret"** might we incur by acting on this evidence?

## 🎯 Key Concept: Deficiency (δ)

The **deficiency** δ measures the information gap between your observational data and a perfect randomized trial:

| δ Value | Interpretation | Action |
|---------|----------------|--------|
| δ ≈ 0 | RCT-quality evidence | Proceed with confidence |
| δ < 0.05 | Excellent | Strong causal conclusions |
| 0.05 ≤ δ < 0.15 | Moderate | Document limitations |
| δ ≥ 0.15 | High | Caution; seek more evidence |

The **safety floor** = 2 × M × δ bounds worst-case policy regret, where M is the utility range.

## Installation

Install the development version from [GitHub](https://github.com/denizakdemir/causaldef):
```r
# install.packages("devtools")
devtools::install_github("denizakdemir/causaldef")
```

## Core Features

| Feature | Function | Description |
|---------|----------|-------------|
| **Deficiency Estimation** | `estimate_deficiency()` | Compare adjustment strategies (IPTW, AIPW, TMLE, etc.) |
| **Policy Regret Bounds** | `policy_regret_bound()` | Compute the "safety floor" for decision-making |
| **Negative Control Diagnostics** | `nc_diagnostic()` | Falsification tests using auxiliary outcomes |
| **Confounding Frontiers** | `confounding_frontier()` | Visualize sensitivity to unmeasured confounding |
| **Survival Analysis** | `causal_spec_survival()` | Time-to-event outcomes (RMST, Cox IPTW) |
| **Front-Door Identification** | `frontdoor_effect()` | When mediators are available |
| **Instrumental Variables** | `iv_effect()` | 2SLS, Wald, and LIML estimators |

## Quick Start: The 4-Step Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│  1. SPECIFY: causal_spec() / causal_spec_survival()             │
│              ↓ Define treatment, outcome, covariates            │
├─────────────────────────────────────────────────────────────────┤
│  2. ESTIMATE: estimate_deficiency()                             │
│               ↓ Compare adjustment methods, select best δ       │
├─────────────────────────────────────────────────────────────────┤
│  3. DIAGNOSE: nc_diagnostic() + confounding_frontier()          │
│               ↓ Test assumptions, map sensitivity               │
├─────────────────────────────────────────────────────────────────┤
│  4. DECIDE: policy_regret_bound() + estimate_effect()           │
│             ↓ Compute safety floor, report causal effect        │
└─────────────────────────────────────────────────────────────────┘
```

## Example 1: Basic Deficiency Estimation

```r
library(causaldef)
set.seed(42)

# Simulate confounded data (W satisfies back-door criterion)
n <- 500
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))  # Treatment confounded by W
Y <- 1 + 2 * A + W + rnorm(n)       # True effect = 2
df <- data.frame(W = W, A = A, Y = Y)

# 1. Define the causal problem
spec <- causal_spec(
  data = df,
  treatment = "A",
  outcome = "Y",
  covariates = "W"
)
#> ✔ Created causal specification: n=500, 1 covariate(s)

# 2. Estimate Le Cam deficiency for different strategies
results <- estimate_deficiency(
  spec, 
  methods = c("unadjusted", "iptw", "aipw"),
  n_boot = 100
)

print(results)
#> ── Le Cam Deficiency Estimates ──────────────────────────────────
#> 
#>      Method  Delta     SE               CI
#>  unadjusted 0.2531 0.0460 [0.1779, 0.3452]
#>        iptw 0.0011 0.0061  [1e-04, 0.0211]
#>        aipw 0.0009 0.0029  [2e-04, 0.0108]
#> 
#> Best method: aipw (δ = 0.0009)
```

**Interpretation:** The unadjusted analysis has δ ≈ 0.25 (substantial confounding). After IPTW or AIPW adjustment, δ drops to ~0.001, indicating near-RCT quality evidence.

## Example 2: Policy Regret Bounds

If we use this evidence to make a policy decision (e.g., approve a drug), what is the worst-case loss?

```r
# Calculate bounds for a utility range of [0, 1]
bounds <- policy_regret_bound(results, utility_range = c(0, 1))
#> ℹ Safety floor: 0.0018 (minimum regret given δ = 0.0009)

print(bounds)
#> ── Policy Regret Bound (Theorem 3.2) ────────────────────────────
#> 
#> • Deficiency δ: 0.0009 
#> • Utility range: [0, 1]
#> • Safety floor: 0.0018 (minimum regret given δ)
#> 
#> Interpretation: Worst-case regret is 0.2% of utility range

plot(bounds, type = "safety_curve")
```

<img src="man/figures/README-regret-1.png" width="100%" />

**The Safety Floor Formula:**

```
Safety Floor = 2 × M × δ
             = 2 × 1 × 0.0009
             = 0.0018
```

This means even with perfect decision-making, the observational evidence introduces at most 0.18% worst-case error.

## Example 3: Negative Control Diagnostic

Test whether adjustment actually removes confounding using a negative control outcome (known to be unaffected by treatment):

```r
# Add a negative control: affected by confounder W, NOT by treatment A
df$Y_nc <- W + rnorm(n)

spec_nc <- causal_spec(
  data = df, 
  treatment = "A", 
  outcome = "Y",
  covariates = "W",
  negative_control = "Y_nc"
)

# Run diagnostic (Theorem 5.2)
nc_test <- nc_diagnostic(spec_nc, method = "iptw")
#> ✔ No evidence against causal assumptions (p = 0.85)

print(nc_test)
#> ── Negative Control Diagnostic ──────────────────────────────────
#> 
#> • δ_NC (observable): 0.0089 
#> • δ bound (Theorem 5.2): 0.0089 (κ = 1)
#> • p-value: 0.85
#> 
#> ✔ NOT FALSIFIED: No evidence against causal assumptions
```

**Decision Logic:**
- p > 0.05 → Adjustment appears adequate
- p ≤ 0.05 → Residual confounding detected; consider additional covariates

## Example 4: Survival Analysis

```r
library(causaldef)
data(hct_outcomes)  # HCT registry data

# Specify survival outcome
spec_surv <- causal_spec_survival(
  data = hct_outcomes,
  treatment = "conditioning_intensity",
  time = "time_to_event",
  event = "event_status",
  covariates = c("age", "disease_status", "kps", "donor_type"),
  estimand = "RMST",
  horizon = 24  # 24-month restricted mean survival
)

# Estimate deficiency
def_surv <- estimate_deficiency(spec_surv, methods = c("unadjusted", "iptw"))

# Compute safety floor in months
bounds_surv <- policy_regret_bound(def_surv, utility_range = c(0, 24))
print(bounds_surv)
#> Safety floor: 1.4 months (worst-case decision error)
```

## Vignettes

Comprehensive tutorials are available:

| Vignette | Description |
|----------|-------------|
| `vignette("introduction")` | Getting started with causaldef |
| `vignette("complete_workflow")` | End-to-end: Specify → Estimate → Diagnose → Decide |
| `vignette("sensitivity_analysis")` | Deficiency vs. E-values comparison |
| `vignette("negative_controls")` | Falsification diagnostics |
| `vignette("survival_analysis")` | Time-to-event outcomes |
| `vignette("policy_learning")` | Safe policy decisions |
| `vignette("causaldef_methodology")` | Theoretical foundations |
| `vignette("advanced_analysis")` | TMLE, matching, GRF methods |

## Theory

Based on Akdemir (2026), ["Constraints on Causal Inference as Experiment Comparison"](https://doi.org/10.5281/zenodo.18367347).

**Core Theorem (Policy Regret Bound):**

```
Regret_do(π) ≤ Regret_obs(π) + 2M × δ
```

Where:
- `Regret_do(π)` = regret under true interventional distribution
- `Regret_obs(π)` = regret estimated from observational data  
- `M` = range of the utility function
- `δ` = Le Cam deficiency

This provides a rigorous justification for using observational evidence in high-stakes decision making, provided δ is small.

## Citation

```bibtex
@misc{causaldef,
  title = {causaldef: Decision-Theoretic Causal Diagnostics via Le Cam Deficiency},
  author = {Akdemir, Deniz},
  year = {2026},
  doi = {10.5281/zenodo.18367347},
  url = {https://github.com/denizakdemir/causaldef},
  note = {R package version 0.2.0}
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

MIT © Deniz Akdemir
