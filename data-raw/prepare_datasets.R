# =============================================================================
# Data Generation Script for causaldef Package
# =============================================================================

set.seed(2025)

# -----------------------------------------------------------------------------
# 1. HCT Outcomes (Survival & Competing Risks)
# -----------------------------------------------------------------------------
n_hct <- 800

# Covariates
age <- rnorm(n_hct, 50, 12)
disease_status <- sample(c("Early", "Intermediate", "Advanced"), n_hct, prob = c(0.4, 0.4, 0.2), replace = TRUE)
kps <- pmin(100, pmax(60, rnorm(n_hct, 85, 10))) # Performance status
donor_type <- sample(c("HLA-Matched", "Mismatched", "Unrelated"), n_hct, replace = TRUE)

# Latent health status (confounder)
U_health <- -0.02 * age + 0.05 * kps - 0.5 * (disease_status == "Advanced") + rnorm(n_hct)

# Treatment: Conditioning Intensity (Myeloablative vs Reduced)
# Reduced (1) often given to older/sicker patients
prob_reduced <- plogis(-2 + 0.05 * age - 0.02 * kps + 0.5 * (disease_status == "Advanced") - 0.5 * U_health)
conditioning <- rbinom(n_hct, 1, prob_reduced)
treatment_label <- ifelse(conditioning == 1, "Reduced", "Myeloablative")

# Outcomes
# Hazard for Death (Event 1)
# Reduced conditioning has slightly higher relapse risk but lower toxicity mortality early on
haz_death <- exp(-4 + 0.3 * conditioning - 0.5 * U_health + 0.02 * age)

# Hazard for Relapse (Event 2 - Competing)
haz_relapse <- exp(-3 + 0.5 * conditioning + 0.3 * (disease_status == "Advanced"))

# Simulate times
time_death <- rexp(n_hct, rate = haz_death)
time_relapse <- rexp(n_hct, rate = haz_relapse)
time_censor <- runif(n_hct, 1, 100) # Administrative censoring

# Observed Data
time_event <- pmin(time_death, time_relapse, time_censor)
event_status <- character(n_hct)
event_status[time_death < time_relapse & time_death < time_censor] <- "Death"
event_status[time_relapse < time_death & time_relapse < time_censor] <- "Relapse"
event_status[event_status == ""] <- "Censored"

hct_outcomes <- data.frame(
  id = 1:n_hct,
  age = round(age),
  disease_status = factor(disease_status),
  kps = round(kps),
  donor_type = factor(donor_type),
  conditioning_intensity = factor(treatment_label),
  time_to_event = round(time_event, 2),
  event_status = factor(event_status)
)

# -----------------------------------------------------------------------------
# 2. Gene Perturbation (Continuous Outcome)
# -----------------------------------------------------------------------------
n_gene <- 500

# Technical factors
batch <- sample(1:4, n_gene, replace = TRUE)
library_size <- rnorm(n_gene, 1e6, 1e5)
cell_cycle <- rnorm(n_gene) # Unobserved confounder
batch_effect <- rnorm(4)[batch]

# Treatment: CRISPR Knockout
# Knockout efficiency affected by library size and cell cycle
prob_ko <- plogis(-1 + 0.5 * scale(library_size) + 0.3 * cell_cycle)
knockout <- rbinom(n_gene, 1, prob_ko)

# Outcomes
# Target Gene: Affected by KO and batch confounders
target_expr <- 10 - 2 * knockout + scale(library_size) + batch_effect + 1.5 * cell_cycle + rnorm(n_gene)

# Housekeeping Gene (Negative Control): NOT affected by KO, but shares confounders
housekeeping_expr <- 8 + scale(library_size) + batch_effect + 1.5 * cell_cycle + rnorm(n_gene)

gene_perturbation <- data.frame(
  sample_id = paste0("S", 1:n_gene),
  batch = factor(batch),
  library_size = round(library_size),
  knockout_status = factor(ifelse(knockout == 1, "Knockout", "Control")),
  target_expression = round(as.numeric(target_expr), 3),
  housekeeping_gene = round(as.numeric(housekeeping_expr), 3)
)

# -----------------------------------------------------------------------------
# Save Data
# -----------------------------------------------------------------------------

# Use usethis to save efficiently if available, else save
if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(hct_outcomes, gene_perturbation, overwrite = TRUE)
} else {
  save(hct_outcomes, file = "data/hct_outcomes.rda", compress = "xz")
  save(gene_perturbation, file = "data/gene_perturbation.rda", compress = "xz")
}
