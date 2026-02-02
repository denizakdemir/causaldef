# =============================================================================
# CausalDef Demo: Complete Workflow Examples
# =============================================================================
# This script demonstrates all major features of the causaldef package.
# Run interactively or source the entire file.
# =============================================================================

library(causaldef)

# Set seed for reproducibility
set.seed(2026)

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  CAUSALDEF: Decision-Theoretic Causal Diagnostics\n")
cat("  Demo Script - All Features\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# =============================================================================
# 1. BASIC DEFICIENCY ESTIMATION
# =============================================================================

cat("\n--- 1. Basic Deficiency Estimation ---\n\n")

# Simulate confounded data
n <- 500
U <- rnorm(n)  # Unmeasured confounder
W <- 0.6 * U + rnorm(n, sd = 0.5)  # Observed proxy
A <- rbinom(n, 1, plogis(0.3 + 0.5 * U))  # Confounded treatment
Y <- 2 + 1.5 * A + 0.8 * U + rnorm(n)  # Outcome

df <- data.frame(W = W, A = A, Y = Y)

# Create causal specification
spec <- causal_spec(
  data = df,
  treatment = "A",
  outcome = "Y",
  covariates = "W"
)

# Estimate deficiency with multiple methods
cat("Estimating deficiency...\n")
def_result <- estimate_deficiency(
  spec,
  methods = c("unadjusted", "iptw", "aipw"),
  n_boot = 100
)

print(def_result)

# Plot deficiency
if (interactive()) {
  print(plot(def_result, type = "bar"))
}

# =============================================================================
# 2. NEGATIVE CONTROL DIAGNOSTIC
# =============================================================================

cat("\n--- 2. Negative Control Diagnostic ---\n\n")

# Add a negative control outcome (affected by U, not by A)
Y_nc <- 0.5 + 0.7 * U + rnorm(n, sd = 0.8)
df$Y_nc <- Y_nc

# Create spec with negative control
spec_nc <- causal_spec(
  data = df,
  treatment = "A",
  outcome = "Y",
  covariates = "W",
  negative_control = "Y_nc"
)

# Run diagnostic
cat("Running negative control diagnostic...\n")
nc_result <- nc_diagnostic(
  spec_nc,
  method = "iptw",
  alpha = 0.05,
  n_boot = 50
)

print(nc_result)

# =============================================================================
# 3. POLICY REGRET BOUNDS
# =============================================================================

cat("\n--- 3. Policy Regret Bounds ---\n\n")

# Calculate safety floor for policy deployment
cat("Computing policy regret bounds...\n")
policy_bounds <- policy_regret_bound(
  deficiency = def_result,
  utility_range = c(0, 10),  # Outcome range
  obs_regret = 0.5  # Observed regret from policy evaluation
)

print(policy_bounds)

# Plot safety curve
if (interactive()) {
  print(plot(policy_bounds, type = "safety_curve"))
}

# =============================================================================
# 4. CONFOUNDING FRONTIER
# =============================================================================

cat("\n--- 4. Confounding Frontier (Sensitivity Analysis) ---\n\n")

# Generate confounding frontier
cat("Computing confounding frontier...\n")
frontier <- confounding_frontier(
  spec,
  alpha_range = c(-2, 2),
  gamma_range = c(-2, 2),
  grid_size = 25
)

print(frontier)

# Percentage of safe region
safe_pct <- mean(frontier$grid$delta < 0.1) * 100
cat(sprintf("\nSafe region (δ < 0.1): %.1f%% of confounding space\n", safe_pct))

# Plot frontier
if (interactive()) {
  print(plot(frontier, type = "heatmap", threshold = c(0.05, 0.1, 0.2)))
}

# =============================================================================
# 5. FRONT-DOOR IDENTIFICATION
# =============================================================================

cat("\n--- 5. Front-Door Identification ---\n\n")

# Simulate front-door scenario
# A -> M -> Y, with U confounding A and Y but not M
M <- 0.5 + 1.2 * A + rnorm(n, sd = 0.5)  # Mediator
Y_fd <- 1 + 0.8 * M + 0.5 * U + rnorm(n)  # Outcome via mediator

df_fd <- data.frame(A = A, M = M, Y = Y_fd)
spec_fd <- causal_spec(df_fd, "A", "Y", covariates = NULL)

# Estimate front-door effect
cat("Estimating front-door effect...\n")
fd_result <- frontdoor_effect(
  spec_fd,
  mediator = "M",
  method = "plugin",
  n_boot = 50
)

print(fd_result)

# =============================================================================
# 6. TRANSPORT DEFICIENCY
# =============================================================================

cat("\n--- 6. Transport Deficiency (Distribution Shift) ---\n\n")

# Source population (e.g., RCT)
n_source <- 400
age_s <- rnorm(n_source, mean = 50, sd = 10)
A_s <- rbinom(n_source, 1, 0.5)  # Randomized
Y_s <- 10 + 2 * A_s - 0.1 * age_s + rnorm(n_source)
source_df <- data.frame(age = age_s, A = A_s, Y = Y_s)

# Target population (different demographics)
n_target <- 200
age_t <- rnorm(n_target, mean = 65, sd = 8)  # Older
target_df <- data.frame(age = age_t)

# Create spec and compute transport
source_spec <- causal_spec(source_df, "A", "Y", "age")

cat("Computing transport deficiency...\n")
transport <- transport_deficiency(
  source_spec,
  target_data = target_df,
  transport_vars = "age"
)

print(transport)

# =============================================================================
# 7. SURVIVAL ANALYSIS
# =============================================================================

cat("\n--- 7. Survival Analysis ---\n\n")

if (requireNamespace("survival", quietly = TRUE)) {
  # Simulate survival data
  n_surv <- 400
  W_s <- rnorm(n_surv)
  A_s <- rbinom(n_surv, 1, plogis(0.3 * W_s))
  
  # Exponential survival times
  lambda <- exp(-0.5 * A_s + 0.2 * W_s)
  time <- rexp(n_surv, rate = lambda)
  censor_time <- rexp(n_surv, rate = 0.3)
  event <- as.integer(time <= censor_time)
  time <- pmin(time, censor_time)
  
  surv_df <- data.frame(W = W_s, A = A_s, time = time, event = event)
  
  # Create survival spec
  surv_spec <- causal_spec_survival(
    surv_df, "A", "time", "event", "W",
    estimand = "RMST",
    horizon = 2
  )
  
  print(surv_spec)
  
  # Estimate deficiency
  cat("\nEstimating survival deficiency...\n")
  surv_def <- estimate_deficiency(
    surv_spec,
    methods = c("iptw", "aipw"),
    n_boot = 50
  )
  
  print(surv_def)
} else {
  cat("Skipping survival analysis (survival package not installed)\n")
}

# =============================================================================
# 8. COMPETING RISKS
# =============================================================================

cat("\n--- 8. Competing Risks ---\n\n")

if (requireNamespace("survival", quietly = TRUE)) {
  # Simulate competing risks: event 1 = death, event 2 = transplant
  n_cr <- 400
  W_cr <- rnorm(n_cr)
  A_cr <- rbinom(n_cr, 1, plogis(0.3 * W_cr))
  
  rate_death <- exp(-0.5 * A_cr + 0.2 * W_cr)
  rate_transplant <- exp(0.3 * A_cr - 0.1 * W_cr)
  
  time_death <- rexp(n_cr, rate_death)
  time_transplant <- rexp(n_cr, rate_transplant)
  time_censor <- runif(n_cr, 0, 3)
  
  obs_time <- pmin(time_death, time_transplant, time_censor)
  event_type <- ifelse(obs_time == time_death, 1,
                       ifelse(obs_time == time_transplant, 2, 0))
  
  cr_df <- data.frame(W = W_cr, A = A_cr, time = obs_time, event = event_type)
  
  # Create competing risks spec
  cr_spec <- causal_spec_competing(
    cr_df, "A", "time", "event", "W",
    event_of_interest = 1,
    horizon = 2
  )
  
  print(cr_spec)
  
  # Estimate deficiency
  cat("\nEstimating competing risks deficiency...\n")
  cr_def <- estimate_deficiency_competing(
    cr_spec,
    method = "cshr",
    n_boot = 30
  )
  
  print(cr_def)
} else {
  cat("Skipping competing risks (survival package not installed)\n")
}

# =============================================================================
# 9. INSTRUMENTAL VARIABLES
# =============================================================================

cat("\n--- 9. Instrumental Variables ---\n\n")

# Simulate IV setting
n_iv <- 800
U_iv <- rnorm(n_iv)  # Unmeasured confounder
Z <- rbinom(n_iv, 1, 0.5)  # Instrument (randomized)
A_iv <- as.numeric(0.3 + 0.4 * Z + 0.3 * U_iv + rnorm(n_iv, sd = 0.3) > 0.5)
Y_iv <- 1 + 2 * A_iv + 0.8 * U_iv + rnorm(n_iv)

iv_df <- data.frame(Z = Z, A = A_iv, Y = Y_iv)
iv_spec <- causal_spec(iv_df, "A", "Y", instrument = "Z")

# Estimate IV effect
cat("Estimating IV effect (LATE)...\n")
iv_result <- iv_effect(
  iv_spec,
  method = "2sls",
  n_boot = 50
)

print(iv_result)

# Test instrument validity
cat("\nTesting instrument validity...\n")
test_instrument(iv_spec)

# =============================================================================
# 10. COMPLETE WORKFLOW SUMMARY
# =============================================================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  WORKFLOW SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Best deficiency achieved: ", round(min(def_result$estimates), 3), "\n")
cat("Best method: ", names(which.min(def_result$estimates)), "\n")
cat("NC diagnostic falsified: ", nc_result$falsified, "\n")
cat("Safety floor (policy): ", round(policy_bounds$safety_floor, 3), "\n")
cat("Transport deficiency: ", round(transport$delta_transport, 3), "\n")

cat("\n")
if (min(def_result$estimates) < 0.1 && !nc_result$falsified) {
  cat("✓ RECOMMENDATION: Causal inference appears reliable.\n")
  cat("  Deploy policy with safety floor = ", round(policy_bounds$safety_floor, 3), "\n")
} else {
  cat("⚠ RECOMMENDATION: Additional caution required.\n")
  cat("  Consider stronger adjustment or randomized trial.\n")
}

cat("\n--- Demo Complete ---\n\n")
