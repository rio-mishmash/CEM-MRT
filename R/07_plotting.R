# ==============================================================================
# R/07_plotting.R: Visualization Functions for Simulation Results
# ==============================================================================

#' Build Custom Comparison Plot (Metrics + Importance)
#' @param sim_output List containing 'Summary' and 'Importance' data frames
#' @param title_str Title for the plot
#' @param xlab_str Label for the X-axis (N (matched))
#' @return A patchwork object combining metrics and importance plots
build_custom_plot <- function(sim_output, title_str = "", xlab_str = "N (matched)") {

  res_df <- sim_output$Summary
  imp_df <- sim_output$Importance
  
  # Prepare metric data
  metric_names <- c("MD_Ratio", "Bias_Ratio", "CI_Ratio", "Coverage", "MSE_Ratio")
  metric_labels <- c("1. MD", "2. Bias", "3. CIW", "4. Coverage", "5. MSE")
  
  plot_sub <- res_df %>%
    filter(Method != "Naive") %>%
    mutate(
      Method = factor(Method, levels = c("PSM", "CEM", "CART", "RF", "MRT")),
      across(all_of(metric_names), as.numeric)
    )
  
  metric_df <- plot_sub %>%
    pivot_longer(
      cols = all_of(metric_names),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      MetricLabel = recode(
        Metric,
        "MD_Ratio" = "1. MD",
        "Bias_Ratio" = "2. Bias",
        "CI_Ratio" = "3. CIW",
        "Coverage" = "4. Coverage",
        "MSE_Ratio" = "5. MSE"
      ),
      MetricLabel = factor(MetricLabel, levels = metric_labels)
    )
  
  x_lims <- if (exists("Nobs", inherits = TRUE)) {
    n_obs <- get("Nobs", inherits = TRUE)
    c(n_obs * 0.60, n_obs)
  } else {
    NULL
  }
  
  # Create metric panels
  p_list <- lapply(seq_along(metric_labels), function(i) {
    metric_i <- metric_labels[i]
    sub_df <- metric_df %>% filter(MetricLabel == metric_i)
    y_lims <- if (metric_i == "3. CIW") c(90, 160) else c(0, 100)
    
    ggplot(sub_df, aes(x = Mean_N, y = Value, color = Method, shape = Method)) +
      geom_line(linewidth = 1.4, alpha = 0.5) +
      geom_point(size = 1.8, alpha = 0.5) +
      facet_wrap(~ MetricLabel) +
      scale_color_manual(
        name = "Method",
        values = c(
          "PSM" = "steelblue",
          "CEM" = "limegreen",
          "CART" = "gold",
          "RF" = "orange",
          "MRT" = "red"
        ),
        guide = guide_legend(nrow = 1)
      ) +
      scale_shape_manual(
        name = "Method",
        values = c("PSM" = 18, "CEM" = 19, "CART" = 18, "RF" = 18, "MRT" = 18),
        guide = "none"
      ) +
      coord_cartesian(xlim = x_lims, ylim = y_lims) +
      labs(
        x = xlab_str,
        y = if (i == 1) "Ratio / Proportion (%)" else ""
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 1)
      )
  })
  
  p_metrics_graphs <- wrap_plots(p_list, ncol = 5)
  
  p_metrics <- wrap_plots(
    p_metrics_graphs,
    guide_area(),
    ncol = 1,
    heights = c(10, 1)
  ) +
    plot_layout(guides = "collect")
  
  # Create variable-importance panel
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    n_all <- unique(imp_df$Target_N)[1]
    
    long_imp <- imp_df %>%
      filter(near(Target_N, n_all * 0.80)) %>%
      pivot_longer(
        cols = starts_with("X"),
        names_to = "Variable",
        values_to = "Value"
      ) %>%
      filter(Method %in% c("CART", "MRT"))
    
    p_imp_graph <- ggplot(
      long_imp,
      aes(x = Variable, y = Value, fill = Method)
    ) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      scale_fill_manual(
        name = "Method",
        values = c("CART" = "gold", "MRT" = "red"),
        guide = guide_legend(nrow = 1)
      ) +
      labs(x = "", y = "Importance") +
      facet_wrap(~ "Importance") +
      coord_flip(ylim = c(0, 1)) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = margin(t = 5, r = 10, b = 5, l = 0)
      )
    
    p_imp <- wrap_plots(
      p_imp_graph,
      guide_area(),
      ncol = 1,
      heights = c(10, 1)
    ) +
      plot_layout(guides = "collect")
  } else {
    p_imp_graph <- ggplot() + theme_void()
    p_imp <- wrap_plots(
      p_imp_graph,
      guide_area(),
      ncol = 1,
      heights = c(10, 1)
    )
  }
  
  # Combine both blocks
  final_plot <- wrap_plots(
    p_metrics,
    p_imp,
    ncol = 2,
    widths = c(5, 1)
  )
  
  return(final_plot)
}



#' Build Model Dependence Violin Plot
build_model_dependence_plot <- function(full_dep_res) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  full_dep_long <- as.data.frame(full_dep_res) %>%
    tidyr::pivot_longer(cols = c("PSM", "MRT"), names_to = "Method", values_to = "Estimate_Range") %>%
    dplyr::filter(!is.na(Estimate_Range))
  
  ggplot(full_dep_long, aes(x = Method, y = Estimate_Range, fill = Method)) +
    geom_violin(alpha = 0.3, trim = TRUE) +
    scale_y_continuous(limits = c(0, 3.6)) +
    facet_wrap(~ Scenario, scales = "free_y") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 12),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title     = element_blank(),
      legend.position  = "right"
    ) +
    labs(y = "Range of estimates", fill = NULL)
}


#' Build Sensitivity Analysis Plot
#' @param res_df Summary data frame from sensitivity simulation
#' @param color_var Name of the parameter to color by (e.g., 'num.trees')
#' @param xlab_str Label for the X-axis
#' @return A patchwork object with sensitivity metrics
build_sensitivity_plot <- function(res_df, color_var, xlab_str = "N (matched)") {
  library(ggplot2)
  library(patchwork)
  
  df1 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$MD_Ratio,   Metric = "1. MD")
  df2 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$Bias_Ratio, Metric = "2. Bias")
  df3 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$CI_Ratio,   Metric = "3. CIW")
  df4 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$Coverage,   Metric = "4. Coverage")
  df5 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$MSE_Ratio,  Metric = "5. MSE")
  long_df <- rbind(df1, df2, df3, df4, df5)
  
  metrics <- unique(long_df$Metric)
  hline_vals <- c("1. MD" = 0, "2. Bias" = 0, "3. CIW" = 100, "4. Coverage" = 100, "5. MSE" = 0)
  p_list <- list()
  
  for (i in 1:5) {
    m <- metrics[i]
    sub_df <- subset(long_df, Metric == m)
    p <- ggplot(sub_df, aes(x = Mean_N, y = Value, color = Color_Val)) +
      geom_point(shape = 16, size = 2.5, alpha = 0.7) +
      facet_wrap(~ Metric) +
      geom_hline(yintercept = hline_vals[m], linetype = "dashed", color = "gray40") +
      scale_color_gradient(low = "red", high = "gold", name = color_var) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title     = element_blank(),
        legend.position = if (i == 5) "right" else "none"
      ) +
      scale_x_continuous(limits = c(160, 200)) +
      labs(y = if(i == 1) "Ratio (%)" else "", x = xlab_str)
    
    if (m == "3. CIW") {
      p <- p + scale_y_continuous(limits = c(100, 150))
    } else {
      p <- p + scale_y_continuous(limits = c(0, 100))
    }
    p_list[[i]] <- p
  }
  
  return(
    wrap_plots(p_list, ncol = 5)
  )
}