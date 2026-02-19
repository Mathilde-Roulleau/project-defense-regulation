# _____Plot clr(ds) vs clr(markers)_____
plot_clr <- function(proteins, ds, markers, markers_clr, condition = NA, use_all = FALSE) {
  df <- correlation_df(proteins, markers_clr, condition, use_all)
  p <- ggplot(df) +
    geom_point(aes(x = median_marker, y = median_ds), color = "grey40", size = 0.5) +
    geom_smooth(aes(x = median_marker, y = median_ds), method = "lm", se = FALSE, color = "grey") +
    labs(
      x = paste0("Median ", markers, " CLR"),
      y = paste0("Median ", ds, " CLR")) +
    annotate("text",
             x = min(df$median_marker),
             y = max(df$median_ds),
             label = condition,
             hjust = 0, vjust = 1,
             size = 3) + 
    theme_minimal()
  
  p
}

plots <- list()

for (cond in conditions) {
  for (ds_ in unique(ds$type)) {
    
    proteins <- ds %>%
      filter(type == ds_) %>%
      pull(protein_in_syst)
    
    # Ribosomal markers
    p_ribo <- plot_clr(
      proteins,
      ds_,
      markers = "ribosomal genes",
      markers_clr = median_growth_markers_clr,
      condition = cond
    )
    
    plots[[paste(cond, ds_, "ribosomal", sep = "_")]] <- p_ribo
    
    # Phage markers
    p_phage <- plot_clr(
      proteins,
      ds_,
      markers = "structural phage genes",
      markers_clr = median_growth_markers_clr,
      condition = cond
    )
    
    plots[[paste(cond, ds_, "phage", sep = "_")]] <- p_phage
  }
}


ncol <- 3
n_plots <- length(plots)
nrow <- ceiling(n_plots / ncol)

combined_plot <- wrap_plots(plots, ncol = ncol)

width_per_col <- 4
height_per_row <- 3

ggsave(
  paste0("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/", project, "/plots_correlations.pdf"),
  combined_plot,
  width  = ncol * width_per_col,
  height = nrow * height_per_row
)