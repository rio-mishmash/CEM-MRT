# ==============================================================================
# R/01_setup.R: Package Management and Parallel Computing Setup
# ==============================================================================

# 1. Package Installation and Loading
# Check for missing packages based on the REQUIRED_PACKAGES list in config.R
new_pkgs <- REQUIRED_PACKAGES[!(REQUIRED_PACKAGES %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  message(paste("Installing missing packages:", paste(new_pkgs, collapse = ", ")))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org/")
}

# Load libraries
suppressPackageStartupMessages({
  # REQUIRED_PACKAGES is defined in config.R
  lapply(REQUIRED_PACKAGES, library, character.only = TRUE)
  library(parallel)
})

# 2. Parallel Backend Configuration
# Detect available physical CPU cores to maximize performance via future.apply
n_cores <- parallel::detectCores(logical = FALSE)
plan(multisession, workers = n_cores)

# 3. Environment & Reproducibility
# Set seed for all subsequent operations
set.seed(GLOBAL_SEED)

# Confirm initialization in the console
message(sprintf(">>> Parallel backend initialized with %d cores.", n_cores))
message(sprintf(">>> ranger (Random Forest) version: %s", packageVersion("ranger")))
message(sprintf(">>> geepack (GEE) version: %s", packageVersion("geepack")))
