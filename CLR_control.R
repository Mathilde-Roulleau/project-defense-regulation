clr_ranked <- function(condition) {
  
  # selection of the columns for the given condition
  cols <- grep(paste0("^", condition, "(_|$)"), colnames(clr), value = TRUE)
  mat <- as.matrix(clr[, cols, drop = FALSE])
  
  # computation of the median and sd on the replicates
  median_clr <- apply(mat, 1, median)
  sd_clr     <- apply(mat, 1, sd)
  
  ribo_flag    <- rownames(clr) %in% rownames(growth_markers_clr)
  phage_flag    <- rownames(clr) %in% rownames(struct_phage_clr)
  defense_flag <- rownames(clr) %in% ds$protein_in_syst
  
  clr_stats <- clr %>%
    mutate(
      median_clr = median_clr,
      sd_clr     = sd_clr,
      rank       = percent_rank(median_clr) * 100,
      ribosome   = ribo_flag,
      phage  = phage_flag,
      defense    = defense_flag
    )
  
  clr_stats$category <- NA
  clr_stats$category[clr_stats$ribosome] <- "Ribosomal genes"
  clr_stats$category[clr_stats$phage]    <- "Phage genes"
  clr_stats$category[clr_stats$defense]  <- "Defense system genes"
  
  
  ggplot(clr_stats, aes(x = rank, y = median_clr)) +
    geom_line(color = "grey20", linewidth = 1) +
    
    geom_errorbar(
      aes(ymin = median_clr - sd_clr, ymax = median_clr + sd_clr),
      width = 0,
      alpha = 0.1,
      color = "grey70"
    ) +
    
    geom_point(
      data = subset(clr_stats, !is.na(category)),
      aes(color = category),
      size = 0.5
    ) +
    
    scale_color_manual(
      values = c(
        "Ribosomal genes" = "red",
        "Phage genes" = "blue4",
        "Defense system genes" = "cyan4"
      ), name = NULL, drop = TRUE) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(size = 9),
      legend.position = c(0.02, 0.98),   
      legend.text = element_text(size = 7),
      legend.justification = c("left", "top"),
      legend.background = element_blank()
    ) +
    
    labs(
      x = "Genes ranked by median expression (%)",
      y = "Median expression level (CLR)",
      title = condition
    )
  
}

plots <- list()

if ("rep" %in% colnames(SRRs)) {
  conditions <- unique(str_remove(colnames(cnt_matrix), "_[^_]*$"))
} else {
  conditions <- unique(colnames(cnt_matrix))
}

for (cond in conditions) {
  plots[[cond]] <- clr_ranked(cond)
}

ncol <- 3
n_plots <- length(plots)
nrow <- ceiling(n_plots / ncol)

combined_plot <- wrap_plots(plots, ncol = ncol)

width_per_col <- 4
height_per_row <- 3

ggsave(
  "EDA/CLR_ranked.png",
  combined_plot,
  width  = ncol * width_per_col,
  height = nrow * height_per_row
)