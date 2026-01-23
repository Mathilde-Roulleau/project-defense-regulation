clr_ranked <- function(condition) {
  
  # selection of the columns for the given condition
  cols <- grep(paste0("^", condition, "(_|$)"), colnames(clr_bacteria), value = TRUE)
  mat <- as.matrix(clr_bacteria[, cols, drop = FALSE])
  
  # computation of the median and sd on the replicates
  median_clr <- apply(mat, 1, median)
  sd_clr     <- apply(mat, 1, sd)
  
  ribo_flag    <- rownames(clr_bacteria) %in% rownames(growth_markers_clr)
  defense_flag <- rownames(clr_bacteria) %in% ds$protein_in_syst
  
  clr_stats <- clr_bacteria %>%
    mutate(
      median_clr = median_clr,
      sd_clr     = sd_clr,
      rank       = percent_rank(median_clr) * 100,
      ribosome   = ribo_flag,
      defense    = defense_flag
    )
  
  ggplot(clr_stats, aes(x = rank, y = median_clr)) +
    geom_line(color = "grey40", linewidth = 1 ) +
    geom_errorbar(
      aes(ymin = median_clr - sd_clr, ymax = median_clr + sd_clr),
      width = 0,
      alpha = 0.3,
      color = "grey40"
    ) +
    geom_point(data = subset(clr_stats, ribosome), color = "salmon", size = 2) +
    geom_point(data = subset(clr_stats, defense), color = "lightskyblue", size = 2) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(), 
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 9)) +
    annotate("text", x = 90, y = 0.2, label = "Ribosomal genes",
             color = "salmon", size = 3, hjust = 1) +
    annotate("text", x = 5, y = 0.2, label = "Defense system genes",
             color = "lightskyblue", size = 3, hjust = 0) +
    labs(
      x = "Genes ranked by median expression (%)",
      y = "Median expression level (CLR)",
      title = condition)
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
