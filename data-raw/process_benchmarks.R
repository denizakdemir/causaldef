
# Script to download and process classical causal inference benchmarks
# 1. Lalonde's NSW (Experimental + Observational Controls)
# 2. Right Heart Catheterization (RHC) - High dimensional confounding

# Load necessary libraries
if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Package 'foreign' is needed to read .dta files.")
}

output_dir <- "../data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# =============================================================================
# 1. Lalonde's NSW Dataset
# =============================================================================
message("Processing Lalonde NSW dataset...")

# URLs from Dehejia's website
url_nsw <- "http://users.nber.org/~rdehejia/data/nsw_dw.dta"
url_cps <- "http://users.nber.org/~rdehejia/data/cps_controls.dta"
url_psid <- "http://users.nber.org/~rdehejia/data/psid_controls.dta"

# Function to safely read DTA
read_dta_url <- function(url) {
  tmp <- tempfile(fileext = ".dta")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  on.exit(unlink(tmp))
  foreign::read.dta(tmp)
}

tryCatch({
  nsw <- read_dta_url(url_nsw)
  cps <- read_dta_url(url_cps)
  psid <- read_dta_url(url_psid)
  
  # Add sample identifiers
  nsw$sample_id <- ifelse(nsw$treat == 1, "nsw_treated", "nsw_control")
  cps$sample_id <- "cps_control"
  psid$sample_id <- "psid_control"
  
  # CPS/PSID from Dehejia don't have 'treat' column (it's implicit 0), add it
  if(!"treat" %in% names(cps)) cps$treat <- 0
  if(!"treat" %in% names(psid)) psid$treat <- 0
  
  # Ensure columns match
  common_cols <- intersect(names(nsw), names(cps))
  
  # Combine
  nsw_benchmark <- rbind(
    nsw[, c(common_cols, "sample_id")],
    cps[, c(common_cols, "sample_id")],
    psid[, c(common_cols, "sample_id")]
  )
  
  # Standardize column names
  # re78 is outcome, treat is treatment
  # Covariates: age, education, black, hispanic, married, nodegree, re74, re75
  
  # Save
  save(nsw_benchmark, file = file.path(output_dir, "nsw_benchmark.rda"), compress = "xz")
  message("Saved nsw_benchmark.rda")
  
}, error = function(e) {
  warning("Failed to download or process Lalonde data: ", e$message)
})

# =============================================================================
# 2. Right Heart Catheterization (RHC)
# =============================================================================
message("Processing RHC dataset...")

# URL from Vanderbilt Biostats (using hbiostat.org)
url_rhc <- "https://hbiostat.org/data/repo/rhc.csv"

tryCatch({
  tmp_rhc <- tempfile(fileext = ".csv")
  # Use curl if possible to handle redirects/https better
  download.file(url_rhc, tmp_rhc, quiet = TRUE)
  
  rhc <- read.csv(tmp_rhc, stringsAsFactors = TRUE)
  
  # Cleanup and subset if too large? 
  # RHC is ~5735 rows, it's fine.
  
  # Key variables:
  # swang1: Treatment (RHC vs No RHC) -> Convert to binary
  # dth30: Outcome (30 day mortality) -> Convert to binary
  # t3d30: Survival time
  
  # Rename for clarity
  # We will keep original names but maybe ensure 'swang1' is clearly factor or 0/1
  
  # Save
  save(rhc, file = file.path(output_dir, "rhc.rda"), compress = "xz")
  message("Saved rhc.rda")
  
}, error = function(e) {
  warning("Failed to download or process RHC data: ", e$message)
  # Fallback: Try alternative URL if primary fails?
  tryCatch({
    url_rhc_alt <- "http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv"
    download.file(url_rhc_alt, tmp_rhc, quiet = TRUE)
    rhc <- read.csv(tmp_rhc, stringsAsFactors = TRUE)
    save(rhc, file = file.path(output_dir, "rhc.rda"), compress = "xz")
    message("Saved rhc.rda (from alt source)")
  }, error = function(e2) {
    warning("Failed alt source for RHC too.")
  })
})

