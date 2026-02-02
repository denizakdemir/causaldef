pkgname <- "causaldef"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('causaldef')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("causal_spec")
### * causal_spec

flush(stderr()); flush(stdout())

### Name: causal_spec
### Title: Create a Causal Problem Specification
### Aliases: causal_spec

### ** Examples

# Create sample data
n <- 200
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))
Y <- 1 + 2 * A + W + rnorm(n)
df <- data.frame(W = W, A = A, Y = Y)

# Create causal specification
spec <- causal_spec(
  data = df,
  treatment = "A",
  outcome = "Y",
  covariates = "W"
)

print(spec)




cleanEx()
nameEx("causal_spec_competing")
### * causal_spec_competing

flush(stderr()); flush(stdout())

### Name: causal_spec_competing
### Title: Causal Specification for Competing Risks
### Aliases: causal_spec_competing

### ** Examples

# Simulate competing risks: death (1) vs transplant (2)
set.seed(42)
n <- 500
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.3 * W))

# Event times
rate_death <- exp(-0.5 * A + 0.2 * W)
rate_transplant <- exp(0.3 * A - 0.1 * W)

time_death <- rexp(n, rate_death)
time_transplant <- rexp(n, rate_transplant)
time_censor <- runif(n, 0, 3)

observed_time <- pmin(time_death, time_transplant, time_censor)
event <- ifelse(observed_time == time_death, 1,
                ifelse(observed_time == time_transplant, 2, 0))

df <- data.frame(W = W, A = A, time = observed_time, event = event)

spec <- causal_spec_competing(
  df, "A", "time", "event", "W",
  event_of_interest = 1,
  horizon = 2
)
print(spec)




cleanEx()
nameEx("causal_spec_survival")
### * causal_spec_survival

flush(stderr()); flush(stdout())

### Name: causal_spec_survival
### Title: Create a Survival Causal Specification
### Aliases: causal_spec_survival

### ** Examples

# Simulate survival data
n <- 200
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))
time <- rexp(n, rate = exp(-0.5 * A + 0.3 * W))
event <- rbinom(n, 1, 0.8)
df <- data.frame(W = W, A = A, time = time, event = event)

spec <- causal_spec_survival(
  data = df,
  treatment = "A",
  time = "time",
  event = "event",
  covariates = "W",
  estimand = "RMST",
  horizon = 5
)




cleanEx()
nameEx("confounding_frontier")
### * confounding_frontier

flush(stderr()); flush(stdout())

### Name: confounding_frontier
### Title: Map the Confounding Frontier
### Aliases: confounding_frontier

### ** Examples

# Basic frontier
frontier <- confounding_frontier(
  alpha_range = c(-2, 2),
  gamma_range = c(-2, 2),
  grid_size = 25
)

# With data-based parameter estimation
df <- data.frame(A = rnorm(100), Y = rnorm(100), W = rnorm(100))
spec <- causal_spec(df, "A", "Y", "W")
frontier <- confounding_frontier(spec, grid_size = 50)




cleanEx()
nameEx("create_plumber_api")
### * create_plumber_api

flush(stderr()); flush(stdout())

### Name: create_plumber_api
### Title: Create Plumber API for CausalDef
### Aliases: create_plumber_api

### ** Examples

## Not run: 
##D # Create API files
##D create_plumber_api("my_api")
##D 
##D # Run the API
##D plumber::plumb("my_api/plumber.R")$run(port = 8080)
## End(Not run)




cleanEx()
nameEx("estimate_deficiency")
### * estimate_deficiency

flush(stderr()); flush(stdout())

### Name: estimate_deficiency
### Title: Estimate Le Cam Deficiency
### Aliases: estimate_deficiency

### ** Examples

# Create sample data
n <- 200
W <- rnorm(n)
A <- rbinom(n, 1, plogis(0.5 * W))
Y <- 1 + 2 * A + W + rnorm(n)
df <- data.frame(W = W, A = A, Y = Y)

spec <- causal_spec(df, "A", "Y", "W")
results <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 50)
print(results)




cleanEx()
nameEx("estimate_effect")
### * estimate_effect

flush(stderr()); flush(stdout())

### Name: estimate_effect
### Title: Estimate Causal Effects from Deficiency Objects
### Aliases: estimate_effect estimate_effect.deficiency

### ** Examples

## Not run: 
##D spec <- causal_spec(df, "A", "Y", "W")
##D def <- estimate_deficiency(spec, methods = "iptw")
##D effect <- estimate_effect(def)
##D print(effect)
## End(Not run)




cleanEx()
nameEx("frontdoor_effect")
### * frontdoor_effect

flush(stderr()); flush(stdout())

### Name: frontdoor_effect
### Title: Front-Door Adjustment Kernel
### Aliases: frontdoor_effect

### ** Examples

# Simulate front-door scenario
n <- 500
U <- rnorm(n)  # Unmeasured confounder
A <- rbinom(n, 1, plogis(0.5 * U))
M <- 0.5 + 1.2 * A + rnorm(n, sd = 0.5)  # Mediator (unconfounded by U)
Y <- 1 + 0.8 * M + 0.5 * U + rnorm(n)    # Outcome
df <- data.frame(A = A, M = M, Y = Y)

spec <- causal_spec(df, "A", "Y", covariates = NULL)
## Not run: 
##D fd_result <- frontdoor_effect(spec, mediator = "M")
##D print(fd_result)
## End(Not run)




cleanEx()
nameEx("iv_effect")
### * iv_effect

flush(stderr()); flush(stdout())

### Name: iv_effect
### Title: Instrumental Variable Effect Estimation
### Aliases: iv_effect

### ** Examples

# Simulate IV setting
set.seed(42)
n <- 1000
U <- rnorm(n)  # Unmeasured confounder
Z <- rbinom(n, 1, 0.5)  # Instrument (randomized encouragement)
A <- 0.3 + 0.4 * Z + 0.3 * U + rnorm(n, sd = 0.3)  # Treatment (continuous)
A <- as.numeric(A > 0.5)  # Dichotomize
Y <- 1 + 2 * A + 0.8 * U + rnorm(n)  # Outcome

df <- data.frame(Z = Z, A = A, Y = Y)
spec <- causal_spec(df, "A", "Y", instrument = "Z")

## Not run: 
##D iv_result <- iv_effect(spec)
##D print(iv_result)
## End(Not run)




cleanEx()
nameEx("nc_diagnostic")
### * nc_diagnostic

flush(stderr()); flush(stdout())

### Name: nc_diagnostic
### Title: Negative Control Diagnostic
### Aliases: nc_diagnostic

### ** Examples

# Create data with negative control
n <- 200
U <- rnorm(n)
W <- U + rnorm(n, sd = 0.5)
A <- rbinom(n, 1, plogis(0.5 * W))  # Binary treatment
Y <- 1 + 2 * A + U + rnorm(n)
Y_nc <- U + rnorm(n)  # Shares U but no effect from A
df <- data.frame(W = W, A = A, Y = Y, Y_nc = Y_nc)

spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")




cleanEx()
nameEx("policy_regret_bound")
### * policy_regret_bound

flush(stderr()); flush(stdout())

### Name: policy_regret_bound
### Title: Compute Policy Regret Bounds
### Aliases: policy_regret_bound

### ** Examples

# From a deficiency estimate
# From a deficiency estimate
df <- data.frame(W=rnorm(100), A=rbinom(100,1,0.5), Y=rnorm(100))
spec <- causal_spec(df, "A", "Y", "W")
def <- estimate_deficiency(spec, methods = "iptw", n_boot = 0)
bound <- policy_regret_bound(def, utility_range = c(0, 1))

# From a numeric value
bound <- policy_regret_bound(
  deficiency = 0.1,
  utility_range = c(0, 100),
  obs_regret = 5
)
print(bound)




cleanEx()
nameEx("run_causaldef_app")
### * run_causaldef_app

flush(stderr()); flush(stdout())

### Name: run_causaldef_app
### Title: Launch CausalDef Shiny Dashboard
### Aliases: run_causaldef_app

### ** Examples

## Not run: 
##D # Launch with example data
##D df <- data.frame(
##D   W = rnorm(200),
##D   A = rbinom(200, 1, 0.5),
##D   Y = rnorm(200)
##D )
##D run_causaldef_app(data = df)
## End(Not run)




cleanEx()
nameEx("transport_deficiency")
### * transport_deficiency

flush(stderr()); flush(stdout())

### Name: transport_deficiency
### Title: Transport Deficiency Between Source and Target Populations
### Aliases: transport_deficiency

### ** Examples

# Source population (RCT)
set.seed(42)
n_source <- 500
age_s <- rnorm(n_source, 50, 10)
A_s <- rbinom(n_source, 1, 0.5)  # Randomized
Y_s <- 10 + 2 * A_s - 0.1 * age_s + rnorm(n_source)
source_df <- data.frame(age = age_s, A = A_s, Y = Y_s, S = 1)

# Target population (different age distribution)
n_target <- 300
age_t <- rnorm(n_target, 65, 8)  # Older population
target_df <- data.frame(age = age_t, S = 0)

source_spec <- causal_spec(source_df, "A", "Y", "age")
## Not run: 
##D transport <- transport_deficiency(
##D   source_spec,
##D   target_data = target_df,
##D   transport_vars = "age"
##D )
##D print(transport)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
