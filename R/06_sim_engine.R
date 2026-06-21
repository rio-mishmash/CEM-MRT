# ==============================================================================
# R/06_sim_engine.R: Simulation Engine and Model Dependence Analysis
# ==============================================================================

#' Run a Simulation Block for a Specific Scenario
#' @param gen_func Data generation function for the scenario (from 02_data_generation.R)
#' @return A list containing Summary data frame and Importance data frame
run_simulation_block <- function(gen_func) {
  # Seed for the simulation block iterations
  set.seed(GLOBAL_SEED)
  
  # 1. Execute Iterations in Parallel using future.apply
  sim_raw <- suppressWarnings(
    future_lapply(1:N_SIM, function(s) {
      tryCatch({
        # Generate dataset for this iteration
        df_sim <- gen_func(N_OBS)
        df_sim$id <- 1:nrow(df_sim)
        x_vars <- grep("^X", names(df_sim), value = TRUE)
        ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars, collapse = " + ")))
        
        # Baseline: Naive (Unadjusted) Estimation
        res_naive <- run_naive(df_sim)
        
        # Baseline: PS Estimation for PSM and standard methods
        m_ps_glm <- glm(formula = ps_formula, data = df_sim, family = binomial(link = "logit"))
        df_sim$ps.est <- predict(m_ps_glm, type = "response")
        
        # Baseline: Standard CEM (fixed 2-bin breaks)
        res_cem <- run_standard_cem(df_sim, breaks_n = 2)
        
        # Pre-train Tree Models for CEM Variants
        # CART (Single Tree)
        base_tree_cart <- ranger::ranger(ps_formula, data = df_sim, splitrule = "gini", num.threads = 1,
                                         min.node.size = MIN_NODE_SIZE, num.trees = 1,
                                         mtry = length(x_vars), max.depth = NULL, node.stats = TRUE)
        
        # Random Forest
        base_rf <- ranger::ranger(ps_formula, data = df_sim, splitrule = "gini", num.threads = 1,
                                  min.node.size = MIN_NODE_SIZE, num.trees = NUM_TREES_RF,
                                  mtry = MTRY_DEFAULT, max.depth = MAX_DEPTH_RF, node.stats = TRUE)
        
        # Forest for MRT Selection
        forest_mrt <- ranger::ranger(ps_formula, data = df_sim, splitrule = "gini", num.threads = 1,
                                     min.node.size = MIN_NODE_SIZE, num.trees = NUM_TREES_MRT,
                                     mtry = MTRY_DEFAULT, max.depth = NULL, node.stats = TRUE)
        
        # Select the Most Representative Tree (MRT)
        rep_tree_obj <- get_representative_tree(forest_mrt, df_sim)
        
        # Results container for current iteration across target sample sizes
        iter_res <- list()
        
        # 2. Iterate through Target Sample Size Grid (defined in config.R)
        for (t_n in TARGET_GRID) {
          # Process CART
          t_info_cart <- prepare_tree_info(base_tree_cart, 1)
          nodes_cart  <- predict(base_tree_cart, df_sim, type = "terminalNodes")$predictions
          res_pruned_cart <- prune_forest_to_target(nodes_cart, t_info_cart, df_sim, t_n)
          imp_cart <- calculate_pruned_importance(base_tree_cart, 1, t_info_cart, res_pruned_cart$final_nodes)
          
          # Process RF
          t_info_rf <- prepare_tree_info(base_rf, 1:NUM_TREES_RF)
          nodes_rf  <- predict(base_rf, df_sim, type = "terminalNodes")$predictions
          res_pruned_rf <- prune_forest_to_target(nodes_rf, t_info_rf, df_sim, t_n)
          imp_rf <- calculate_pruned_importance(base_rf, 1:NUM_TREES_RF, t_info_rf, res_pruned_rf$final_nodes)
          
          # Process MRT
          t_info_mrt <- prepare_tree_info(forest_mrt, rep_tree_obj$index)
          nodes_mrt  <- predict(forest_mrt, df_sim, type = "terminalNodes")$predictions[, rep_tree_obj$index, drop = FALSE]
          res_pruned_mrt <- prune_forest_to_target(nodes_mrt, t_info_mrt, df_sim, t_n)
          imp_mrt <- calculate_pruned_importance(forest_mrt, rep_tree_obj$index, t_info_mrt, res_pruned_mrt$final_nodes)
          
          # Store all method results for this Target N
          iter_res[[as.character(t_n)]] <- list(
            Naive = res_naive,
            PSM   = run_psm_subclass_target(df_sim, t_n, already_estimated = TRUE),
            CART  = calc_gee_stratum_ate(dplyr::mutate(df_sim, s = res_pruned_cart$strata), "s"),
            RF    = calc_gee_stratum_ate(dplyr::mutate(df_sim, s = res_pruned_rf$strata), "s"),
            MRT   = calc_gee_stratum_ate(dplyr::mutate(df_sim, s = res_pruned_mrt$strata), "s"),
            CEM   = res_cem,
            Imp_CART = imp_cart,
            Imp_RF   = imp_rf,
            Imp_MRT  = imp_mrt
          )
        }
        return(iter_res)
      }, error = function(e) { return(NULL) })
    }, future.seed = TRUE)
  )
  
  # 3. Aggregate Raw Results into Summary Tables
  return(aggregate_sim_results(sim_raw))
}

#' Helper to Extract Metrics Matrix for a Specific Target and Method
get_metrics_mat <- function(raw_list, target, method) {
  extracted <- lapply(raw_list, function(x) x[[as.character(target)]][[method]])
  # Filter out NULL results from failed iterations
  do.call(rbind, extracted[!sapply(extracted, is.null)])
}

#' Helper to Aggregate Simulation Results and Calculate Metrics Ratios
aggregate_sim_results <- function(sim_raw) {
  res_summary_list <- list()
  imp_raw_list <- list()
  idx <- 1
  imp_idx <- 1
  
  # Retrieve baseline reference values from the Naive method at the first grid point
  # Fix: Pass only the first element [1] instead of the entire TARGET_GRID
  naive_mat_base <- get_metrics_mat(sim_raw, TARGET_GRID[1], "Naive")
  naive_md_base   <- mean(as.numeric(naive_mat_base[, "md"]), na.rm = TRUE)
  naive_bias_base <- mean(as.numeric(naive_mat_base[, "bias"]), na.rm = TRUE)
  naive_mse_base  <- mean(as.numeric(naive_mat_base[, "bias"])^2, na.rm = TRUE)
  naive_ci_base   <- mean(as.numeric(naive_mat_base[, "ci_width"]), na.rm = TRUE)
  
  for (m in METHODS) {
    for (t_n in TARGET_GRID) {
      metrics <- get_metrics_mat(sim_raw, t_n, m)
      if (is.null(metrics) || nrow(metrics) == 0) next
      
      # Calculate Relative Ratios (%) compared to Naive performance
      res_summary_list[[idx]] <- data.frame(
        Method     = m,
        Target_N   = t_n,
        Mean_N     = mean(as.numeric(metrics[, "n_valid"]), na.rm = TRUE),
        MD_Ratio   = (mean(as.numeric(metrics[, "md"]), na.rm = TRUE) / (if(naive_md_base == 0) 1 else naive_md_base)) * 100,
        Bias_Ratio = (abs(mean(as.numeric(metrics[, "bias"]), na.rm = TRUE)) / (if(abs(naive_bias_base) == 0) 1 else abs(naive_bias_base))) * 100,
        CI_Ratio   = (mean(as.numeric(metrics[, "ci_width"]), na.rm = TRUE) / (if(naive_ci_base == 0) 1 else naive_ci_base)) * 100,
        MSE_Ratio  = (mean(as.numeric(metrics[, "bias"])^2, na.rm = TRUE) / (if(naive_mse_base == 0) 1 else naive_mse_base)) * 100
      )
      idx <- idx + 1
      
      # Aggregate Variable Importance for Tree methods
      if (m %in% c("CART", "RF", "MRT")) {
        imp_list <- lapply(sim_raw, function(x) x[[as.character(t_n)]][[paste0("Imp_", m)]])
        imp_df   <- as.data.frame(do.call(rbind, imp_list[!sapply(imp_list, is.null)]))
        if (nrow(imp_df) > 0) {
          imp_df$Method <- m
          imp_df$Target_N <- t_n
          imp_raw_list[[imp_idx]] <- imp_df
          imp_idx <- imp_idx + 1
        }
      }
    }
  }
  
  return(list(
    Summary    = do.call(rbind, res_summary_list),
    Importance = if (length(imp_raw_list) > 0) do.call(rbind, imp_raw_list) else NULL
  ))
}

#' Run Model Dependence Simulation
#' Evaluates the robustness of estimates across various covariate adjustment sets
#' @param gen_func Data generation function for the specific scenario
#' @param n_obs Number of observations for data generation
#' @param n_sim Number of simulation iterations
#' @param target_n Target matched sample size for stratification
#' @return A data frame containing the range (max - min) of ATE estimates for PSM and MRT
run_model_dependence_simulation <- function(gen_func, n_obs, n_sim, target_n) {
  # 1. Define all possible covariate adjustment sets (2^5 = 32 patterns)
  all_covs <- c("X1", "X2", "X3", "X4", "X5")
  test_patterns <- list()
  for (k in 1:length(all_covs)) {
    combos <- combn(all_covs, k, simplify = FALSE)
    for (comb in combos) {
      name <- paste(comb, collapse = "+")
      test_patterns[[name]] <- comb
    }
  }
  # Add a baseline pattern with no additional covariate adjustment in the GEE
  test_patterns[["None"]] <- NULL
  
  # 2. Execute parallel simulation loop
  dep_results <- future_lapply(1:n_sim, function(s) {
    tryCatch({
      # Generate dataset for the current iteration
      df <- gen_func(n_obs)
      df$id <- 1:nrow(df)
      
      # Propensity Score estimation remains constant (using all base variables)
      x_vars_all <- grep("^X", names(df), value = TRUE)
      ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars_all, collapse = " + ")))
      
      # Fit PS model and calculate scores
      m_ps_glm <- glm(formula = ps_formula, data = df, family = binomial(link = "logit"))
      df$ps.est <- predict(m_ps_glm, type = "response")
      
      # Train forest and select the Most Representative Tree (MRT)
      forest_mrt <- ranger::ranger(ps_formula, data = df, splitrule = "gini", num.threads = 1,
                                   num.trees = 100, mtry = 3, min.node.size = 0.05 * nrow(df), 
                                   replace = TRUE)
      rep_tree_obj <- get_representative_tree(forest_mrt, df)
      rep_tree_obj$model <- forest_mrt
      
      # 3. Estimate ATE for every covariate pattern using PSM (Subclassification)
      psm_ests <- sapply(names(test_patterns), function(name) {
        res <- run_psm_subclass_target(df, target_n, already_estimated = TRUE, 
                                       include_covariates = test_patterns[[name]])
        return(as.numeric(res["est"]))
      })
      
      # 4. Estimate ATE for every covariate pattern using MRT (CEM-based)
      mrt_ests <- sapply(names(test_patterns), function(name) {
        res <- run_cem_mrt_target(df, target_n, rep_tree = rep_tree_obj, 
                                  include_covariates = test_patterns[[name]])
        return(as.numeric(res["est"]))
      })
      
      # 5. Calculate Model Dependence as the range (Max - Min) of estimates
      p_range <- if(all(is.na(psm_ests))) NA else max(psm_ests, na.rm = TRUE) - min(psm_ests, na.rm = TRUE)
      m_range <- if(all(is.na(mrt_ests))) NA else max(mrt_ests, na.rm = TRUE) - min(mrt_ests, na.rm = TRUE)
      
      data.frame(Iter = s, PSM = p_range, MRT = m_range)
      
    }, error = function(e) return(NULL))
  }, future.seed = TRUE)
  
  # Return combined results, excluding failed iterations
  return(do.call(rbind, dep_results[!sapply(dep_results, is.null)]))
}

#' Run MRT Hyperparameter Sensitivity Analysis
#' Evaluates how changing num.trees, mtry, and min.node.size affects performance
#' @param gen_func Data generation function for the scenario
#' @param fixed_t_n The specific matched sample size to evaluate (e.g., 180)
#' @param param_grid A data frame grid of hyperparameters to test
#' @return A data frame containing metric ratios relative to a naive baseline
run_sensitivity_simulation_block <- function(gen_func, fixed_t_n, param_grid = NULL) {
  set.seed(GLOBAL_SEED)
  
  # 1. Execute parallel simulation iterations
  sim_raw <- suppressWarnings(
    future_lapply(1:N_SIM, function(s) {
      tryCatch({
        # Generate data for the current iteration
        df_sim <- gen_func(N_OBS)
        df_sim$id <- 1:nrow(df_sim)
        
        # Calculate naive baseline for this specific dataset
        res_naive <- run_naive(df_sim)
        
        iter_res <- list()
        
        # 2. Iterate through the hyperparameter grid
        for (i in 1:nrow(param_grid)) {
          p_ntree <- param_grid$num.trees[i]
          p_mtry <- param_grid$mtry[i]
          p_minsize <- param_grid$min.node.size[i]
          
          # Create a unique key for this configuration
          p_key <- paste0("T:", p_ntree, ",M:", p_mtry, ",N:", p_minsize)
          
          # Run MRT CEM for the current hyperparameter set
          iter_res[[p_key]] <- list(
            Naive = res_naive, 
            MRT = run_cem_mrt_target(
              df_sim, 
              fixed_t_n, 
              num.trees = p_ntree, 
              mtry = p_mtry, 
              min.node.size = p_minsize
            )
          )
        }
        return(iter_res)
      }, error = function(e) { NULL })
    }, future.seed = TRUE)
  )
  
  # Filter out failed iterations
  sim_raw <- sim_raw[!sapply(sim_raw, is.null)]
  
  # Internal helper to extract metrics from the raw simulation list
  get_metrics_mat <- function(raw_list, p_key, method) {
    extracted <- lapply(raw_list, function(x) x[[p_key]][[method]])
    do.call(rbind, extracted[!sapply(extracted, is.null)])
  }
  
  # 3. Calculate Baseline Reference Values (from Naive estimates)
  # We use the first key in the grid to extract the shared naive results
  first_key <- paste0("T:", param_grid$num.trees[1], ",M:", param_grid$mtry[1], ",N:", param_grid$min.node.size[1])
  naive_mat <- get_metrics_mat(sim_raw, first_key, "Naive")
  
  naive_md_base   <- mean(as.numeric(naive_mat[, "md"]), na.rm = TRUE)
  naive_bias_base <- mean(as.numeric(naive_mat[, "bias"]), na.rm = TRUE)
  naive_mse_base  <- mean(as.numeric(naive_mat[, "bias"])^2, na.rm = TRUE)
  naive_ci_base   <- mean(as.numeric(naive_mat[, "ci_width"]), na.rm = TRUE)
  
  # 4. Aggregate Results and Calculate Ratios (%)
  res_list <- list()
  for (i in 1:nrow(param_grid)) {
    p_ntree   <- param_grid$num.trees[i]
    p_mtry    <- param_grid$mtry[i]
    p_minsize <- param_grid$min.node.size[i]
    p_key     <- paste0("T:", p_ntree, ",M:", p_mtry, ",N:", p_minsize)
    
    metrics <- get_metrics_mat(sim_raw, p_key, "MRT")
    
    res_list[[i]] <- data.frame(
      Method        = p_key, 
      num.trees     = p_ntree, 
      mtry          = p_mtry, 
      min.node.size = p_minsize, 
      Target_N      = fixed_t_n,
      Mean_N        = mean(as.numeric(metrics[, "n_valid"]), na.rm = TRUE),
      MD_Ratio      = (mean(as.numeric(metrics[, "md"]), na.rm = TRUE) / (if(naive_md_base == 0) 1 else naive_md_base)) * 100,
      Bias_Ratio    = (abs(mean(as.numeric(metrics[, "bias"]), na.rm = TRUE)) / (if(abs(naive_bias_base) == 0) 1 else abs(naive_bias_base))) * 100,
      CI_Ratio      = (mean(as.numeric(metrics[, "ci_width"]), na.rm = TRUE) / (if(naive_ci_base == 0) 1 else naive_ci_base)) * 100,
      MSE_Ratio     = (mean(as.numeric(metrics[, "bias"])^2, na.rm = TRUE) / (if(naive_mse_base == 0) 1 else naive_mse_base)) * 100
    )
  }
  
  return(do.call(rbind, res_list))
}
