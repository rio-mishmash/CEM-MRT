# ==============================================================================
# config.R: Global Parameters and Simulation Settings
# ==============================================================================

# --- General R Environment Options ---
# Disable scientific notation for better readability of results
options(scipen = 999)
# Expand print width and max output
options(width = 1000, max.print = 9999)

# --- Reproducibility ---
# Seed for simulation reproducibility
GLOBAL_SEED <- 2026

# --- Simulation Parameters ---
# Number of observations (sample size) per iteration
N_OBS <- 200 

# Number of simulation iterations
N_SIM <- 1000

# --- Core Method Parameters ---

# PSM Subclass Settings
MAX_Q_BINS <- floor(N_OBS / 4) # Maximum number of subclasses for PSM

# Tree-based CEM Settings (CART, Random Forest, MRT)
TREE_SEED      <- 2026
MIN_NODE_SIZE  <- 0.05 * N_OBS # Minimum node size as 5% of N_OBS
NUM_TREES_RF   <- 10           # Number of trees for Random Forest
NUM_TREES_MRT  <- 100          # Number of trees for Most Representative Tree
MTRY_DEFAULT   <- 2            # Default number of variables to sample for splits
MAX_DEPTH_RF   <- 3            # Tree depth limit for RF

# --- Output Directories ---
# Paths for saving results and visualizations
OUTPUT_DIR     <- "output/"
TABLE_PATH     <- "output/tables/"
PLOT_PATH      <- "output/plots/"

# --- Methodology Constants ---
# List of methods to be evaluated in the simulation engine
METHODS <- c("Naive", "PSM", "CEM", "CART", "RF", "MRT")

# Labels for plotting metrics
METRIC_LABELS <- c("1. MD", "2. Bias", "3. CI Width", "4. MSE")

# --- Package List ---
# Required packages for the simulation
REQUIRED_PACKAGES <- c(
  "ggh4x",        # graph panels
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