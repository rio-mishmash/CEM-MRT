# ==============================================================================
# R/05_methods.R: Wrapper Functions for Causal Inference Methods
# ==============================================================================

#' Naive (Unmatched) Estimation
#' @param df Data frame with Y, Z, and covariates
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results for the full sample
run_naive <- function(df, include_covariates = FALSE) {
  calc_gee_stratum_ate(df, strata_col = NULL, include_covariates = include_covariates)
}

#' Propensity Score Matching (Subclassification) to Target Sample Size
#' @param df Data frame with Y, Z, and covariates
#' @param target_n Target number of matched units
#' @param already_estimated If TRUE, uses existing 'ps.est' column in df
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results after PS subclassification
run_psm_subclass_target <- function(df, target_n, already_estimated = FALSE, include_covariates = FALSE) {
  if (!already_estimated) {
    x_vars <- grep("^X", names(df), value = TRUE)
    ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars, collapse = " + ")))
    m_ps_glm <- glm(formula = ps_formula, data = df, family = binomial(link = "logit"))
    df$ps.est <- predict(m_ps_glm, type = "response")
  }
  
  ps_min <- min(df$ps.est, na.rm = TRUE)
  ps_max <- max(df$ps.est, na.rm = TRUE)
  
  if (!is.finite(ps_min) || !is.finite(ps_max) || ps_min == ps_max) {
    return(c(est = NA, bias = NA, ci_width = NA, coverage = NA, md = NA, max_smd = NA, n_valid = 0, n_strata = 0))
  }
  
  best_q <- 1
  min_diff <- Inf
  max_Q <- floor(nrow(df) / 4)
  
  # Search for optimal number of subclasses to match target_n
  for (Q in max_Q:1) {
    q_cuts <- seq(ps_min, ps_max, length.out = Q + 1)
    temp_strata <- cut(df$ps.est, breaks = q_cuts, include.lowest = TRUE, labels = FALSE)
    tab <- table(factor(temp_strata, levels = 1:Q), df$Z)
    
    if (!"0" %in% colnames(tab) || !"1" %in% colnames(tab)) next
    
    valid_bins <- tab[, "1"] >= 2 & tab[, "0"] >= 2
    current_n  <- sum(rowSums(tab)[valid_bins])
    diff_val   <- abs(current_n - target_n)
    
    if (diff_val < min_diff) {
      min_diff <- diff_val
      best_q <- Q
    }
    if (diff_val <= 5) break
  }
  
  final_cuts <- seq(ps_min, ps_max, length.out = best_q + 1)
  df$ps_strata <- cut(df$ps.est, breaks = final_cuts, include.lowest = TRUE, labels = FALSE)
  
  calc_gee_stratum_ate(df, "ps_strata", include_covariates = include_covariates)
}

#' Standard Coarsened Exact Matching (CEM)
#' @param df Data frame with Y, Z, and covariates
#' @param breaks_n Number of bins for each covariate
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results for CEM strata
run_standard_cem <- function(df, breaks_n, include_covariates = FALSE) {
  df_cem <- df
  # Identify all X variables
  x_vars <- grep("^X", names(df), value = TRUE)
  
  # Create bins for each covariate
  strata_matrix <- matrix(NA, nrow = nrow(df), ncol = length(x_vars))
  for (i in seq_along(x_vars)) {
    strata_matrix[, i] <- cut(df[[x_vars[i]]], breaks = breaks_n, labels = FALSE, include.lowest = TRUE)
  }
  
  # Combine bins into unique strata
  df_cem$cem_strata <- as.integer(factor(apply(strata_matrix, 1, paste, collapse = "_")))
  
  calc_gee_stratum_ate(df_cem, "cem_strata", include_covariates = include_covariates)
}

#' CEM using CART (Single Tree)
#' @param df Data frame with Y, Z, and covariates
#' @param target_n Target number of matched units
#' @param base_tree Optional pre-trained ranger tree (1-tree forest)
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results for CART-based strata
run_cem_cart_target <- function(df, target_n, base_tree = NULL, include_covariates = FALSE) {
  if (is.null(base_tree)) {
    x_vars <- grep("^X", names(df), value = TRUE)
    ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars, collapse = " + ")))
    # Train a single tree using ranger (analogous to CART)
    base_tree <- ranger::ranger(ps_formula, data = df, splitrule = "gini", num.threads = 1,
                                num.trees = 1, min.node.size = 0.05 * nrow(df), 
                                replace = FALSE, sample.fraction = 1, seed = TREE_SEED)
  }
  
  tree_info <- prepare_tree_info(base_tree, 1)
  initial_nodes <- predict(base_tree, df, type = "terminalNodes")$predictions
  if (is.vector(initial_nodes)) initial_nodes <- matrix(initial_nodes, ncol = 1)
  
  # Prune tree to reach target sample size
  res_pruned <- prune_forest_to_target(initial_nodes, tree_info, df, target_n)
  df$cart_strata <- res_pruned$strata
  
  calc_gee_stratum_ate(df, "cart_strata", include_covariates = include_covariates)
}

#' CEM using Random Forest (RF)
#' @param df Data frame with Y, Z, and covariates
#' @param target_n Target number of matched units
#' @param base_rf Optional pre-trained ranger forest
#' @param num.trees Number of trees in the forest
#' @param mtry Number of variables for splitting
#' @param max.depth Maximum depth of trees
#' @param min.node.size Minimum size of terminal nodes
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results for RF-based strata
run_cem_rf_target <- function(df, target_n, base_rf = NULL, 
                              num.trees = 10, mtry = 3, max.depth = 3, 
                              min.node.size = 0.05 * nrow(df), include_covariates = FALSE) {
  if (is.null(base_rf)) {
    x_vars <- grep("^X", names(df), value = TRUE)
    ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars, collapse = " + ")))
    base_rf <- ranger::ranger(ps_formula, data = df, splitrule = "gini", num.threads = 1,
                              num.trees = num.trees, mtry = mtry, max.depth = max.depth, 
                              min.node.size = min.node.size, replace = TRUE, seed = TREE_SEED)
  }
  
  tree_info_list <- prepare_tree_info(base_rf, 1:base_rf$num.trees)
  initial_nodes <- predict(base_rf, df, type = "terminalNodes")$predictions
  
  # Prune forest to reach target sample size
  res_pruned <- prune_forest_to_target(initial_nodes, tree_info_list, df, target_n)
  df$rf_strata <- res_pruned$strata
  
  calc_gee_stratum_ate(df, "rf_strata", include_covariates = include_covariates)
}

#' CEM using Most Representative Tree (MRT)
#' @param df Data frame with Y, Z, and covariates
#' @param target_n Target number of matched units
#' @param rep_tree Optional pre-selected MRT object
#' @param num.trees Number of trees to train forest
#' @param mtry Number of variables for splitting
#' @param min.node.size Minimum size of terminal nodes
#' @param include_covariates Logical or vector of covariate names for GEE adjustment
#' @return GEE results for MRT-based strata
run_cem_mrt_target <- function(df, target_n, rep_tree = NULL, 
                               num.trees = 100, mtry = 3, 
                               min.node.size = 0.05 * nrow(df), include_covariates = FALSE) {
  if (is.null(rep_tree)) {
    x_vars <- grep("^X", names(df), value = TRUE)
    ps_formula <- as.formula(paste("as.factor(Z) ~", paste(x_vars, collapse = " + ")))
    forest <- ranger::ranger(ps_formula, data = df, splitrule = "gini", num.threads = 1,
                             num.trees = num.trees, mtry = mtry,
                             max.depth = NULL, min.node.size = min.node.size, 
                             replace = TRUE, seed = TREE_SEED)
    # Select MRT from the forest
    rep_tree <- get_representative_tree(forest, df)
    rep_tree$model <- forest
  }
  
  tree_info <- prepare_tree_info(rep_tree$model, rep_tree$index)
  # Extract terminal nodes specifically for the MRT
  initial_nodes <- predict(rep_tree$model, df, type = "terminalNodes")$predictions[, rep_tree$index, drop = FALSE]
  
  # Prune MRT to reach target sample size
  res_pruned <- prune_forest_to_target(initial_nodes, tree_info, df, target_n)
  df$mrt_strata <- res_pruned$strata
  
  calc_gee_stratum_ate(df, "mrt_strata", include_covariates = include_covariates)
}
