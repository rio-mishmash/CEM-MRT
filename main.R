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

# Define Scenarios to evaluate
scenarios <- list(
  list(id = "Linear_High_Corr", func = gen_scen_high_corr, title = "Linear - High Correlation"),
  list(id = "Linear_Low_Corr",  func = gen_scen_low_corr,  title = "Linear - Low Correlation"),
  list(id = "Tree_Based",       func = gen_scen_tree,      title = "Tree-based Scenario"),
  list(id = "Complex_Terms",    func = gen_scen_complex,   title = "Complex (Higher-order) Terms"),
  list(id = "HTE",              func = gen_scen_hte,       title = "Heterogeneous Treatment Effect")
)

# 4. Simulation Loop across Sample Sizes (N=200 and N=400)
#for (n_val in c(N_OBS, N_OBS*2)) {
for (n_val in c(N_OBS)) {
  
  message(sprintf("\n>>> Starting Simulation Batch for Sample Size N = %d <<<", n_val))
  
  for (scen in scenarios) {
    # message(sprintf("Running Scenario: %s", scen$title))
    tictoc::tic(sprintf("Scenario: %s (N=%d)", scen$id, n_val))
    
    # Execute simulation engine
    # Note: run_simulation_block uses future_lapply for parallel execution
    res <- run_simulation_block(scen$func, 
                                run_methods = c("PSM", "CEM", "CART"), 
                                cutoff_percentiles=NULL) #seq(0.1, 0.9, by=0.2)
    
    if (!is.null(res$Summary)) {
      # Store results in the global list
      result_key <- paste0(scen$id, "_N", n_val)
      all_results[[result_key]] <- res$Summary
      
      # Create output directories if needed
      if (!dir.exists(TABLE_PATH)) dir.create(TABLE_PATH, recursive = TRUE)
      if (!dir.exists(PLOT_PATH))  dir.create(PLOT_PATH,  recursive = TRUE)
      
      # Define output file paths
      csv_file  <- paste0(TABLE_PATH, "summary_", result_key, ".csv")
      plot_file <- paste0(PLOT_PATH,  "plot_",    result_key, ".png")
      
      # Remove existing files
      if (file.exists(csv_file))  file.remove(csv_file)
      if (file.exists(plot_file)) file.remove(plot_file)
      
      # Export summary table to CSV
      write.csv(res$Summary, file = csv_file, row.names = FALSE)
      
      # Generate and save visualization
      p <- build_custom_plot(res, scen$title)
      ggsave(filename = plot_file, plot = p, width = 12, height = 4)
      
    } else {
      warning(sprintf("Simulation failed for Scenario: %s", scen$id))
    }
    
    tictoc::toc()
  }
}

beepr::beep(sound = 2, expr = NULL)
message("\nWorkflow Complete.")