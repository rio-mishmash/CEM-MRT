# ==============================================================================
# config.R: Global Parameters and Simulation Settings
# ==============================================================================

# --- 1. General R Environment Options ---
# Disable scientific notation for better readability of results
options(scipen = 999)
# Expand print width and max output
options(width = 1000, max.print = 9999)

# --- 2. Reproducibility ---
# Seed for simulation reproducibility
GLOBAL_SEED <- 2026

# --- 3. Simulation Parameters ---
# Number of observations (sample size) per iteration
# Common values: 200 or 400
N_OBS <- 200 

# Number of simulation iterations
N_SIM <- 2000

# --- 4. Target Sample Size Grid ---
# Defines the range of matched sample sizes to evaluate
# Default: From N_OBS down to 50% of N_OBS in 5% decrements
TARGET_GRID <- seq(N_OBS, floor(N_OBS * 0.70), by = -floor(N_OBS * 0.10))

# --- 5. Core Method Parameters ---

# PSM Subclass Settings
MAX_Q_BINS <- floor(N_OBS / 4) # Maximum number of subclasses for PSM

# Tree-based CEM Settings (CART, Random Forest, MRT)
TREE_SEED      <- 2026
MIN_NODE_SIZE  <- 0.05 * N_OBS # Minimum node size as 5% of N_OBS
NUM_TREES_RF   <- 10           # Number of trees for Random Forest
NUM_TREES_MRT  <- 100          # Number of trees for Most Representative Tree
MTRY_DEFAULT   <- 2            # Default number of variables to sample for splits
MAX_DEPTH_RF   <- 3            # Tree depth limit for RF

# --- 6. Output Directories ---
# Paths for saving results and visualizations
OUTPUT_DIR     <- "output/"
TABLE_PATH     <- "output/tables/"
PLOT_PATH      <- "output/plots/"

# --- 7. Methodology Constants ---
# List of methods to be evaluated in the simulation engine
METHODS <- c("Naive", "PSM", "CEM", "CART", "RF", "MRT")

# Labels for plotting metrics
METRIC_LABELS <- c("1. MD", "2. Bias", "3. CI Width", "4. MSE")

# --- 8. Package List ---
# Required packages for the simulation
REQUIRED_PACKAGES <- c(
  "geepack",      # Generalized Estimating Equations (GEE)
  "ranger",       # Fast Random Forest implementation
  "future.apply", # Parallel processing
  "ggplot2",      # Data visualization
  "dplyr",        # Data manipulation
  "MASS",         # Multivariate normal distribution
  "patchwork",    # Plot composition
  "tictoc",       # Execution timing
  "tidyr"         # Data tidying
)