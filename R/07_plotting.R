# ==============================================================================
# R/07_plotting.R: Visualization Functions for Simulation Results
# ==============================================================================

#' Build Custom Comparison Plot (Metrics + Importance)
#' @param sim_output List containing 'Summary' and 'Importance' data frames
#' @param title_str Title for the plot
#' @param xlab_str Label for the X-axis (Matched Sample Size)
#' @return A patchwork object combining metrics and importance plots
build_custom_plot <- function(sim_output, title_str = "", xlab_str = "Matched Sample Size") {
  library(patchwork)
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  
  res_df <- sim_output$Summary
  imp_df <- sim_output$Importance
  
  # 1. Prepare Data for Metric Plots
  plot_sub <- res_df %>% dplyr::filter(Method != "Naive")
  plot_sub$Method <- factor(plot_sub$Method, levels = c("PSM", "CEM", "CART", "RF", "MRT"))
  
  # Metrics as defined in project configuration
  metrics_labels <- c("1. MD", "2. Bias", "3. CI Width", "4. MSE")
  p_list <- list()
  
  for (i in 1:4) {
    val_col <- switch(i, "MD_Ratio", "Bias_Ratio", "CI_Ratio", "MSE_Ratio")
    
    temp_df <- plot_sub
    temp_df$Value <- as.numeric(temp_df[[val_col]])
    temp_df$MetricLabel <- metrics_labels[i]
    
    p <- ggplot(temp_df, aes(x = Mean_N, y = Value, color = Method, shape = Method)) +
      geom_line(linewidth = 1.2, alpha = 0.4) + 
      geom_point(size = 2.4, alpha = 0.6) +
      facet_wrap(~ MetricLabel) + 
      theme_bw() +
      scale_color_manual(
        name = "Method",
        values = c("PSM" = "steelblue", "CEM" = "limegreen", "CART" = "gold", "RF" = "orange", "MRT" = "red"),
        guide = guide_legend(nrow = 1)
      ) +
      scale_shape_manual(
        name = "Method",
        values = c("PSM" = 18, "CEM" = 19, "CART" = 18, "RF" = 18, "MRT" = 18), 
        guide = "none"
      ) +
      labs(y = if(i == 1) "Ratio (%)" else "", x = xlab_str) +
      theme(
        plot.title      = element_text(size = 12, face = "bold"),
        axis.title.x    = element_text(size = 11),
        axis.title.y    = element_text(size = 11),
        axis.text.x     = element_text(size = 11),
        axis.text.y     = element_text(size = 12),
        strip.text      = element_text(size = 12),
        legend.position = "bottom" 
      )
    
    x_lims <- if (exists("Nobs")) c(Nobs * 0.60, Nobs) else NULL
    y_lims <- if (i == 3) c(90, 160) else c(0, 100)
    
    p <- p + coord_cartesian(xlim = x_lims, ylim = y_lims)
    p_list[[i]] <- p
  }
  
  # 2. Prepare Variable Importance Plot
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    target_1 <- unique(imp_df$Target_N)[1]
    long_imp <- imp_df %>% 
      dplyr::filter(Target_N == target_1) %>% 
      tidyr::pivot_longer(cols = starts_with("X"), names_to = "Variable", values_to = "Value") %>% 
      dplyr::filter(Method %in% c("CART", "MRT"))
    
    # Box Plot
    p_imp_graph <- ggplot(long_imp, aes(x = Variable, y = Value, fill = Method)) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      theme_bw() +
      scale_fill_manual(
        name = "Method",
        values = c("CART" = "gold", "RF" = "orange", "MRT" = "red"),
        guide = guide_legend(nrow = 1)
      ) +
      labs(y = "Importance", x = "") +
      facet_wrap(~ "5. Importance") +
      coord_flip(ylim = c(0, 1)) +
      theme(
        axis.text.y     = element_text(size = 12),
        strip.text      = element_text(size = 12),
        legend.position = "bottom"
      )
    
    # ## Violin Plot
    # p_imp_graph <- ggplot(long_imp, aes(x = Variable, y = Value, fill = Method)) +
    #   geom_violin(alpha = 0.7, scale = "area") + 
    #   theme_bw() +
    #   scale_fill_manual(
    #     name = "Method",
    #     values = c("CART" = "gold", "RF" = "orange", "MRT" = "red"), 
    #     guide = guide_legend(nrow = 1)
    #   ) +
    #   labs(y = "Importance", x = "") + 
    #   facet_wrap(~ "5. Importance") +
    #   coord_flip(ylim = c(0, 1)) +
    #   theme(
    #     axis.text.y     = element_text(size = 12),
    #     strip.text      = element_text(size = 12),
    #     legend.position = "bottom" 
    #   )
    
    # Match the exact vertical structure of the metrics block
    p_imp <- patchwork::wrap_plots(
      p_imp_graph,
      patchwork::guide_area(),
      ncol = 1,
      heights = c(10, 1)
    ) + patchwork::plot_layout(guides = "collect")
    
  } else {
    p_imp_graph <- ggplot() + theme_void()
    p_imp <- patchwork::wrap_plots(p_imp_graph, patchwork::guide_area(), ncol = 1, heights = c(10, 1))
  }
  
  # 3. Combine with patchwork (Perfectly synchronized layout)
  
  # Block A: 4 metric plots with a forced guide_area at the bottom
  p_metrics_graphs <- patchwork::wrap_plots(p_list[1:4], ncol = 4)
  
  p_metrics <- patchwork::wrap_plots(
    p_metrics_graphs, 
    patchwork::guide_area(), 
    ncol = 1, 
    heights = c(10, 1)
  ) + 
    patchwork::plot_layout(guides = "collect")
  
  # Block B: The importance plot (p_imp) now has the identical 10:1 structure as Block A
  
  # Final assembly: Because both blocks share the exact same structural dimensions, 
  # their panel heights and widths will align flawlessly when combined horizontally.
  final_plot <- patchwork::wrap_plots(p_metrics, p_imp, ncol = 2, widths = c(4, 1))
    #patchwork::plot_annotation(title = title_str)
  
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
    scale_y_continuous(limits = c(0, 3)) +
    facet_wrap(~ Scenario, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x      = element_text(size = 11),
      axis.text.y      = element_text(size = 12),
      strip.text       = element_text(size = 12),
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
build_sensitivity_plot <- function(res_df, color_var, xlab_str = "Matched Sample Size") {
  library(ggplot2)
  library(patchwork)
  
  df1 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$MD_Ratio,   Metric = "1. MD")
  df2 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$Bias_Ratio, Metric = "2. Bias")
  df3 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$CI_Ratio,   Metric = "3. CI Width")
  df4 <- data.frame(Mean_N = res_df$Mean_N, Color_Val = res_df[[color_var]], Value = res_df$MSE_Ratio,  Metric = "4. MSE")
  long_df <- rbind(df1, df2, df3, df4)
  
  metrics <- unique(long_df$Metric)
  hline_vals <- c("1. MD" = 0, "2. Bias" = 0, "3. CI Width" = 100, "4. MSE" = 0)
  p_list <- list()
  
  for (i in 1:4) {
    m <- metrics[i]
    sub_df <- subset(long_df, Metric == m)
    p <- ggplot(sub_df, aes(x = Mean_N, y = Value, color = Color_Val)) +
      geom_point(shape = 16, size = 2.5, alpha = 0.7) +
      facet_wrap(~ Metric) +
      geom_hline(yintercept = hline_vals[m], linetype = "dashed", color = "gray40") +
      scale_color_gradient(low = "gold", high = "red", name = color_var) +
      theme_bw() +
      theme(
        strip.text      = element_text(size = 11),
        axis.text.y     = element_text(size = 12),
        legend.title     = element_blank(),
        legend.position = "right" # Legend on the right
      ) +
      scale_x_continuous(limits = c(160, 200)) +
      labs(y = if(i == 1) "Ratio (%)" else "", x = xlab_str)
    
    if (m == "3. CI Width") {
      p <- p + scale_y_continuous(limits = c(100, 150))
    } else {
      p <- p + scale_y_continuous(limits = c(0, 100))
    }
    p_list[[i]] <- p
  }
  
  wrap_plots(p_list, ncol = 4) + 
    plot_layout(guides = "collect") + 
    theme(legend.position = "right")
}