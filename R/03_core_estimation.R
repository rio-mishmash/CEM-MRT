# ==============================================================================
# R/03_core_estimation.R: Core Functions for Effect Estimation and Balance
# ==============================================================================

#' Calculate Stratum-Specific ATE using GEE
#' @param df Data frame containing Y, Z, and covariates
#' @param strata_col Name of the column containing strata IDs (NULL for naive)
#' @param include_covariates Logical or character vector of covariates to adjust for in GEE
#' @return A named vector of metrics: est, bias, ci_width, coverage, md, max_smd, n_valid, n_strata
calc_gee_stratum_ate <- function(df, strata_col = NULL, include_covariates = FALSE) {
  
  # 1. Stratum Validation and Weighting
  if (is.null(strata_col)) {
    valid_strata <- 1
    df_valid <- df
    df_valid$weight <- 1
    df_valid$id <- 1:nrow(df_valid)
  } else {
    strata <- df[[strata_col]]
    # Ensure each stratum has at least 2 treated and 2 control units for stable estimation
    tab_t <- tapply(df$Z, strata, function(z) sum(z == 1))
    tab_c <- tapply(df$Z, strata, function(z) sum(z == 0))
    valid_strata <- names(tab_t)[which(tab_t >= 2 & tab_c >= 2)]
    
    df_valid <- df[strata %in% valid_strata, ]
    
    if (nrow(df_valid) == 0) {
      return(c(est = NA, bias = NA, ci_width = NA, coverage = NA, 
               md = NA, max_smd = NA, n_valid = 0, n_strata = 0))
    }
    
    # Calculate Matching Weights
    s_valid <- df_valid[[strata_col]]
    N_s  <- tapply(df_valid$Z, s_valid, length)
    N_ts <- tapply(df_valid$Z, s_valid, sum)
    N_cs <- N_s - N_ts
    
    p_ts <- N_ts[as.character(s_valid)] / N_s[as.character(s_valid)]
    p_cs <- N_cs[as.character(s_valid)] / N_s[as.character(s_valid)]
    
    # Matching Weight formula: min(ps, 1-ps) / ps_actual
    df_valid$weight <- ifelse(df_valid$Z == 1,  pmin(p_ts, p_cs) / p_ts, pmin(p_ts, p_cs) / p_cs)
    
    if (length(valid_strata) >= 2) {
      df_valid$id <- df_valid[[strata_col]]
    } else {
      df_valid$id <- 1:nrow(df_valid)
    }
  }
  
  # 2. Balance Metrics (MD and SMD)
  all_x_vars <- grep("^X", names(df_valid), value = TRUE)
  X_valid <- as.matrix(df_valid[, all_x_vars])
  Z_valid <- df_valid$Z
  W_valid <- df_valid$weight
  
  # Weighted Means
  X_T_w_mean <- colSums(X_valid[Z_valid == 1, , drop = FALSE] * W_valid[Z_valid == 1]) / sum(W_valid[Z_valid == 1])
  X_C_w_mean <- colSums(X_valid[Z_valid == 0, , drop = FALSE] * W_valid[Z_valid == 0]) / sum(W_valid[Z_valid == 0])
  diff_mean  <- X_T_w_mean - X_C_w_mean
  
  # Standardized Mean Difference (SMD)
  sd_full <- apply(as.matrix(df[, all_x_vars]), 2, sd, na.rm = TRUE)
  max_smd <- max(abs(diff_mean / sd_full), na.rm = TRUE)
  
  # Mahalanobis Distance (MD)
  Sigma <- cov(as.matrix(df[, all_x_vars]), use = "pairwise.complete.obs")
  md <- tryCatch({
    as.numeric(sqrt(t(diff_mean) %*% solve(Sigma) %*% diff_mean))
  }, error = function(e) { NA })
  
  # 3. GEE Estimation
  # Determine formula
  if (is.character(include_covariates)) {
    gee_formula <- as.formula(paste("Y ~ Z +", paste(include_covariates, collapse = " + ")))
  } else if (isTRUE(include_covariates)) {
    gee_formula <- as.formula(paste("Y ~ Z +", paste(all_x_vars, collapse = " + ")))
  } else {
    gee_formula <- as.formula("Y ~ Z")
  }
  
  # Fit GEE model
  fit <- tryCatch({
    geepack::geeglm(gee_formula, data = df_valid, weights = weight, id = id, 
                    corstr = "independence", family = gaussian)
  }, error = function(e) { NULL })
  
  if (is.null(fit)) {
    return(c(est = NA, bias = NA, ci_width = NA, coverage = NA, md = md, 
             max_smd = max_smd, n_valid = nrow(df_valid), n_strata = length(valid_strata)))
  }
  
  # Extract results
  true_ate <- if("TE" %in% colnames(df_valid)) mean(df_valid$TE) else NA
  res_est  <- as.numeric(coef(fit)["Z"])
  res_se   <- as.numeric(sqrt(diag(vcov(fit)))["Z"])
  ci_half  <- 1.96 * res_se
  
  # Coverage: 1 if true ATE is within CI
  cov_val <- if (!is.na(true_ate)) {
    as.integer((res_est - ci_half <= true_ate) & (true_ate <= res_est + ci_half))
  } else {
    NA
  }
  
  out <- c(
    est      = res_est,
    bias     = if(!is.na(true_ate)) res_est - true_ate else NA,
    ci_width = 2 * ci_half,
    coverage = cov_val,
    md       = md,
    max_smd  = max_smd,
    n_valid  = nrow(df_valid),
    n_strata = length(valid_strata)
  )
  
  return(round(out, digits = 5))
}