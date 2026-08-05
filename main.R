# ==============================================================================
# main.R: Execution script for the Causal Inference Simulation
# ==============================================================================

# 1. Load configuration and initialize environment
source("config.R")
source("R/01_setup.R")

# 2. Load modular functions
source("R/02_data_generation.R")
source("R/03_core_estimation.R")
source("R/04_tree_utils.R")
source("R/05_methods.R")
source("R/06_sim_engine.R")
source("R/07_plotting.R")

# 3. Initialization
all_results <- list()

# Seed for the simulation block iterations
set.seed(GLOBAL_SEED)

# Define Scenarios to evaluate
scenarios <- list(
  list(id = "Linear_High_Corr", func = gen_scen_high_corr, title = "Linear - High Correlation"),
  list(id = "Linear_Low_Corr",  func = gen_scen_low_corr,  title = "Linear - Low Correlation"),
  list(id = "Tree_Based",       func = gen_scen_tree,      title = "Tree-based Scenario"),
  list(id = "Complex_Terms",    func = gen_scen_complex,   title = "Complex (Higher-order) Terms"),
  list(id = "HTE",              func = gen_scen_hte,       title = "Heterogeneous Treatment Effect")
)

# Create output directories if needed
if (!dir.exists(TABLE_PATH)) dir.create(TABLE_PATH, recursive = TRUE)
if (!dir.exists(PLOT_PATH))  dir.create(PLOT_PATH,  recursive = TRUE)

# Delete old files
unlink(file.path(TABLE_PATH, "*"), recursive = TRUE)
unlink(file.path(PLOT_PATH,  "*"), recursive = TRUE)

# 4. Simulation Loop across Sample Sizes
for (n_val in c(N_OBS, N_OBS*2)) {
  
  message(sprintf("\n>>> Starting Simulation Batch for Sample Size N = %d <<<", n_val))
  
  # Target Sample Size Grid
  TARGET_GRID <- seq(n_val, floor(n_val * 0.60), by = -floor(n_val * 0.05))
  
  for (scen in scenarios) {
    message(sprintf("Running Scenario: %s", scen$title))
    tictoc::tic(sprintf("Scenario: %s (N=%d)", scen$id, n_val))
    
    # Execute simulation engine
    # Note: run_simulation_block uses future_lapply for parallel execution
    res <- run_simulation_block(scen$func, n_val,
                                run_methods = c("PSM", "CEM", "CART", "RF", "MRT"),
                                cutoff_percentiles = seq(0.1, 0.9, by=0.05))
    
    if (!is.null(res$Summary)) {
      # Store results in the global list
      result_key <- paste0(scen$id, "_N", n_val)
      all_results[[result_key]] <- res$Summary
      
      # Export summary table to CSV
      write.csv(res$Summary, 
                file = paste0(TABLE_PATH, "summary_", result_key, ".csv"), 
                row.names = FALSE)
      
      # Generate and save the visualization (Metrics + Importance)
      p <- build_custom_plot(res, scen$title)
      ggsave(filename = paste0(PLOT_PATH, "plot_", result_key, ".png"), 
             plot = p, width = 12, height = 3.6)
      
    } else {
      warning(sprintf("Simulation failed for Scenario: %s", scen$id))
    }
    
    tictoc::toc()
  }
}


# 2. Model Dependence Analysis
message("\n>>> Running Model Dependence Analysis <<<")
res_dep_high <- run_model_dependence_simulation(gen_scen_high_corr, N_OBS, N_SIM, N_OBS*0.8)
res_dep_low  <- run_model_dependence_simulation(gen_scen_low_corr,  N_OBS, N_SIM, N_OBS*0.8)
res_dep_comp <- run_model_dependence_simulation(gen_scen_complex,   N_OBS, N_SIM, N_OBS*0.8)

res_dep_high$Scenario <- "1.High Corr"; 
res_dep_low$Scenario  <- "2.Low Corr"; 
res_dep_comp$Scenario <- "3.Complex"
full_dep_res <- rbind(res_dep_high, res_dep_low, res_dep_comp)

ggsave(paste0(PLOT_PATH, "plot_model_dependence.png"),
       build_model_dependence_plot(full_dep_res), width = 11, height = 3)


# 3. MRT Sensitivity Analysis
message("\n>>> Running MRT Sensitivity Analysis <<<")
gen_scen_sens <- function(n) { # Nonlinear scenario
  df <- generate_base_df(n, rho = 0.2); df$TE <- 3
  df$RS <- 1.0 * (1 * exp(df$X1 - 1) + 1 * ((df$X2 + 1)^2) + 3 * df$X3 + 2 * df$X4 + 1 * df$X5 + 2 * df$X4 * df$X5)
  df$PS <- 1 / (1 + exp(-0.1 * (3 * df$X1 + 1 * df$X2 + 3 * df$X3 + 2 * df$X4 + 1 * df$X5)))
  df$Z <- rbinom(n, 1, df$PS); df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}

# Execute across grids of num.trees, mtry, and min.node.size
res_trees <- run_sensitivity_simulation_block(gen_scen_sens, N_OBS*0.9,
                                              expand.grid(num.trees = c(20, 30, 40, 50, 100), mtry = 3, min.node.size = 20))
ggsave(paste0(PLOT_PATH, "sens_mrt_trees.png"), build_sensitivity_plot(res_trees, "num.trees"), width = 11, height = 3)

res_mtry <- run_sensitivity_simulation_block(gen_scen_sens, N_OBS*0.9,
                                             expand.grid(num.trees = 50, mtry = c(2, 3, 4, 5), min.node.size = 20))
ggsave(paste0(PLOT_PATH, "sens_mrt_mtry.png"), build_sensitivity_plot(res_mtry, "mtry"), width = 11, height = 3)

res_node <- run_sensitivity_simulation_block(gen_scen_sens, N_OBS*0.9,
                                             expand.grid(num.trees = 50, mtry = 3, min.node.size = c(5, 10, 15, 20)))
ggsave(paste0(PLOT_PATH, "sens_mrt_nodesize.png"), build_sensitivity_plot(res_node, "min.node.size"), width = 11, height = 3)


beepr::beep(sound = 3, expr = NULL)
message("\nWorkflow Complete.")