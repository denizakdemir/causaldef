
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causaldef

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/causaldef)](https://CRAN.R-project.org/package=causaldef)
<!-- badges: end -->

**causaldef** implements Le Cam deficiency theory for causal inference,
providing quantitative bounds on information loss from confounding,
selection bias, and distributional shift.

Unlike traditional sensitivity analysis which focuses on “how much bias”
exists, `causaldef` answers the decision-theoretic question: **“how much
regret”** might we incur by acting on this evidence?

## Key Concept: Deficiency (δ)

The **deficiency** δ is a theoretical measure of the information gap
between your observational data and a perfect randomized trial. In
practice, `causaldef` provides a **computable proxy** $\widehat{\delta}$
based on propensity-score TV balance (PS-TV), which is informative about
overlap/positivity and residual confounding risk.

For bounded utilities with range $M$ (max minus min), the manuscript
provides:

- a **regret transfer penalty** (upper bound term) of $M\cdot \delta$,
  and
- a minimax **safety floor** (lower bound) of $(M/2)\cdot \delta$.

`policy_regret_bound()` reports both quantities.

## Installation

You can install the development version of causaldef from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("denizakdemir/causaldef")
```

## Core Features

- **Deficiency proxies:** `estimate_deficiency()` (PS-TV overlap/balance
  proxy)
- **Policy regret bounds:** `policy_regret_bound()` (transfer penalty +
  minimax floor)
- **Negative control diagnostics:** `nc_diagnostic()` (falsification and
  bounds)
- **Sensitivity analysis:** `confounding_frontier()` (linear-Gaussian
  confounding frontier)
- **Survival + competing risks:** `causal_spec_survival()`,
  `causal_spec_competing()`

## Example 1: Basic Deficiency Estimation

``` r
library(causaldef)
set.seed(42)

# Simulate confounded data (W satisfies back-door criterion)
n <- 500
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))
Y <- 1 + 2 * A + W + rnorm(n)
df <- data.frame(W = W, A = A, Y = Y)

# 1. Define the causal problem
spec <- causal_spec(
  data = df,
  treatment = "A",
  outcome = "Y",
  covariates = "W"
)
#> ✔ Created causal specification: n=500, 1 covariate(s)

# 2. Estimate a deficiency proxy (PS-TV) for different strategies
results <- estimate_deficiency(
  spec, 
  methods = c("unadjusted", "iptw", "aipw"),
  n_boot = 100
)
#> ℹ Estimating deficiency: unadjusted
#> ℹ Estimating deficiency: iptw
#> ℹ Estimating deficiency: aipw

print(results)
#> 
#> -- Deficiency Proxy Estimates (PS-TV) ------
#> 
#>      Method  Delta     SE               CI            Quality
#>  unadjusted 0.1190 0.0230 [0.1048, 0.1982] Insufficient (Red)
#>        iptw 0.0212 0.0099 [0.0142, 0.0537]  Excellent (Green)
#>        aipw 0.0212 0.0087 [0.0154, 0.0483]  Excellent (Green)
#> Note: delta is a propensity-score TV proxy (overlap/balance diagnostic).
#> 
#> Best method: aipw (delta = 0.0212 )
```

**Interpretation:** Unadjusted $\widehat{\delta} \approx$ 0.119; after
IPTW/AIPW, $\widehat{\delta} \approx$ 0.021.

## Example 2: Policy Regret Bounds

If we use this evidence to make a policy decision (e.g., approve a
drug), what is the worst-case loss?

``` r
# Calculate bounds for a utility range of [0, 1]
bounds <- policy_regret_bound(results, utility_range = c(0, 1))
#> ℹ Transfer penalty: 0.0212 (delta = 0.0212)

print(bounds)
#> 
#> -- Policy Regret Bounds -------------------------------------------------
#> 
#> * Deficiency delta: 0.0212 
#> * Delta mode: point 
#> * Delta method: aipw 
#> * Utility range: [0, 1]
#> * Transfer penalty: 0.0212 (additive regret upper bound)
#> * Minimax floor: 0.0106 (worst-case lower bound)
#> 
#> Interpretation: Transfer penalty is 2.1 % of utility range given delta
plot(bounds, type = "safety_curve")
```

<img src="man/figures/README-regret-1.png" width="100%" />

The plug-in transfer penalty is 0.0212 on a 0–1 utility scale; the
minimax safety floor is 0.0106.

## Example 3: Negative Control Diagnostic

Check if the “Adjusted” strategy actually removes confounding using a
negative control outcome $Y_{nc}$ (known to be unaffected by treatment).

``` r
# Add a negative control to simulation
df$Y_nc <- W + rnorm(n) # Correlated with W (confounder) but not A

spec_nc <- causal_spec(
  data = df, 
  treatment = "A", 
  outcome = "Y",
  covariates = "W",
  negative_control = "Y_nc"
)
#> ✔ Created causal specification: n=500, 1 covariate(s)

# Run diagnostic
nc_test <- nc_diagnostic(spec_nc, method = "iptw")
#> ℹ Using kappa = 1 (conservative). Consider domain-specific estimation or sensitivity analysis via kappa_range.
#> ✔ No evidence against causal assumptions (p = 0.85072 )
print(nc_test)
#> 
#> -- Negative Control Diagnostic ----------------------------------------
#> 
#> * delta_NC (observable): 0.0089 
#> * delta bound (NC bound): 0.0089 (kappa = 1 )
#> * p-value: 0.85072 
#> 
#> RESULT: NOT REJECTED. We failed to catch an error, but that doesn't mean an error isn't there.
#> NOTE: Your effect estimate must exceed the Noise Floor (delta_bound) to be meaningful.
```

Here, the test does not reject (p = 0.851), and the observable proxy is
$\widehat{\delta}_{NC} \approx$ 0.009.

## Example 4: Survival Analysis (HCT)

``` r
data(hct_outcomes)

# Create an explicit 0/1 event indicator (any non-censor event)
hct <- hct_outcomes
hct$event_any <- as.integer(hct$event_status != "Censored")

spec_surv <- causal_spec_survival(
  data = hct,
  treatment = "conditioning_intensity",
  time = "time_to_event",
  event = "event_any",
  covariates = c("age", "disease_status", "kps", "donor_type"),
  estimand = "RMST",
  horizon = 24
)
#> ✔ Created survival causal specification: n=800, 677 events

def_surv <- estimate_deficiency(spec_surv, methods = c("unadjusted", "cox_iptw"), n_boot = 50)
#> ℹ Inferred treatment value: Reduced
#> ℹ Estimating deficiency: unadjusted
#> ℹ Estimating deficiency: cox_iptw
print(def_surv)
#> 
#> -- Deficiency Proxy Estimates (PS-TV) ------
#> 
#>      Method  Delta     SE               CI            Quality
#>  unadjusted 0.3030 0.0596 [0.2254, 0.4282] Insufficient (Red)
#>    cox_iptw 0.0076 0.0047  [0.0079, 0.024]  Excellent (Green)
#> Note: delta is a propensity-score TV proxy (overlap/balance diagnostic).
#> 
#> Best method: cox_iptw (delta = 0.0076 )

bounds_surv <- policy_regret_bound(def_surv, utility_range = c(0, 24))
#> ℹ Transfer penalty: 0.1823 (delta = 0.0076)
print(bounds_surv)
#> 
#> -- Policy Regret Bounds -------------------------------------------------
#> 
#> * Deficiency delta: 0.0076 
#> * Delta mode: point 
#> * Delta method: cox_iptw 
#> * Utility range: [0, 24]
#> * Transfer penalty: 0.1823 (additive regret upper bound)
#> * Minimax floor: 0.0911 (worst-case lower bound)
#> 
#> Interpretation: Transfer penalty is 0.8 % of utility range given delta
```

## Theory

Based on Akdemir (2026), [“Constraints on Causal Inference as Experiment
Comparison”](https://doi.org/10.5281/zenodo.18367347).

The core theorem links the deficiency $\delta$ (Total Variation
distance) to the max-min regret:

$$ \text{Regret}_{do}(\pi) \leq \text{Regret}_{obs}(\pi) + M \cdot \delta $$

Where $M$ is the range of the utility function. In practice, the package
provides plug-in bounds by feeding a computable proxy/estimate (e.g.,
$\widehat{\delta}$) into the regret formula.

## Citation

``` bibtex
@misc{causaldef,
  title = {causaldef: Decision-Theoretic Causal Diagnostics via Le Cam Deficiency},
  author = {Akdemir, Deniz},
  year = {2026},
  doi = {10.5281/zenodo.18367347}
}
```
