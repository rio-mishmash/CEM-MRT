# ==============================================================================
# R/02_data_generation.R: Functions for Simulation Scenarios
# ==============================================================================

#' Generate Base Covariates
#' @param n Sample size
#' @param rho Correlation coefficient among variables
#' @return A data frame with 5 covariates (X1-X2 continuous, X3-X5 binary)
generate_base_df <- function(n, rho = 0) {
  mu <- rep(0, 5)
  Sigma <- matrix(rho, nrow = 5, ncol = 5)
  diag(Sigma) <- 1
  
  # Generate multivariate normal distribution
  X <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)
  
  # Convert to data frame and round/binarize
  df <- data.frame(
    X1 = round(X[, 1], digits = 2),
    X2 = round(X[, 2], digits = 2),
    X3 = round(X[, 3], digits = 2),
    X4 = round(X[, 4], digits = 2),
    X5 = round(X[, 5], digits = 2)
    # Transform continuous variables to binary (-1 or 1) based on sign
    # X1 = sign(X[, 1]),
    # X2 = sign(X[, 2]),
    # X3 = sign(X[, 3]),
    # X4 = sign(X[, 4]),
    # X5 = sign(X[, 5])
  )
  return(df)
}

#' Scenario: Linear - High Correlation
gen_scen_high_corr <- function(n) {
  df <- generate_base_df(n, rho = 0.2)
  df$TE <- 3
  df$RS    <- 1.0 * ( 3 * df$X1 + 2 * df$X2 + 1 * df$X3 + 2 * df$X4 + 1 * df$X5)
  logit_ps <- 0.1 * ( 6 * df$X1 + 4 * df$X2 + 2 * df$X3 - 1 * df$X4 + 1 * df$X5)
  df$PS <- 1 / (1 + exp(-logit_ps))
  df$Z <- rbinom(n, 1, df$PS)
  df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}

#' Scenario: Linear - Low Correlation
gen_scen_low_corr <- function(n) {
  df <- generate_base_df(n, rho = 0.2)
  df$TE <- 3
  df$RS    <- 1.0 * ( 3 * df$X1 + 2 * df$X2 + 1 * df$X3 + 2 * df$X4 + 1 * df$X5)
  logit_ps <- 0.1 * ( 6 * df$X1 - 4 * df$X2 + 2 * df$X3 - 1 * df$X4 + 1 * df$X5)
  df$PS <- 1 / (1 + exp(-logit_ps))
  df$Z <- rbinom(n, 1, df$PS)
  df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}

#' Scenario: Tree-Based Propensity Score
gen_scen_tree <- function(n) {
  df <- generate_base_df(n, rho = 0.2)
  df$TE <- 3
  df$RS    <- 1.0 * ( 3 * df$X1 + 2 * df$X2 + 1 * df$X3 + 2 * df$X4 + 1 * df$X5)
  # Step-function (tree-like) propensity score logic
  df$PS <- dplyr::case_when(
    df$X3 > 0 & df$X4 > 0 ~ 0.80,
    df$X3 > 0 & df$X4 < 0 & df$X1 >  -1 ~ 0.60,
                            df$X1 <= -1 ~ 0.20,
    TRUE ~ 0.40
  )
  df$Z <- rbinom(n, 1, df$PS)
  df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}

#' Scenario: Complex (Higher-order) Terms
gen_scen_complex <- function(n) {
  df <- generate_base_df(n, rho = 0.2)
  df$TE <- 3
  df$RS    <- 1.0 * ( 3 * df$X1 + 2 * df$X2 + 1 * df$X3 + 2 * df$X4 + 1 * df$X5)
  logit_ps <- 0.1 * ( 6 * df$X1 + 4 * df$X2 + 2 * df$X3 - 1 * df$X4 + 1 * df$X5  +
                      1 * ((df$X1 + 1)^2) - exp(df$X2 - 1) - 2 * df$X4 * df$X5 - 1)
  df$PS <- 1 / (1 + exp(-logit_ps))
  df$Z <- rbinom(n, 1, df$PS)
  df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}

#' Scenario: Heterogeneous Treatment Effect (HTE)
gen_scen_hte <- function(n) {
  df <- generate_base_df(n, rho = 0.2)
  # Treatment Effect varies based on covariates
  df$TE    <- 3 + 
              0.1 * ( 3 * df$X1 - 2 * df$X2 + 1 * df$X3 - 2 * df$X4 + 1 * df$X5)
  df$RS    <- 1.0 * ( 3 * df$X1 + 2 * df$X2 + 1 * df$X3 + 2 * df$X4 + 1 * df$X5)
  logit_ps <- 0.1 * ( 6 * df$X1 + 4 * df$X2 + 2 * df$X3 - 1 * df$X4 + 1 * df$X5)
  df$PS <- 1 / (1 + exp(-logit_ps))
  df$Z <- rbinom(n, 1, df$PS)
  df$Y <- df$RS + df$TE * df$Z + rnorm(n, 0, 1)
  return(df)
}
