## Compute PCoA (based on CLR values) and plot it

dist_matrix_clr <- vegdist(t(clr), method = "euclidean")

# Compute principal coordinate decomposition
pcoa_clr <- pcoa(dist_matrix_clr)

pcoa_df_clr <- data.frame(
  condition = rownames(pcoa_clr$vectors),
  Axis1 = pcoa_clr$vectors[, 1],
  Axis2 = pcoa_clr$vectors[, 2]
)

pcoa_df_clr <- left_join(pcoa_df_clr, SRRs, by = "condition")

# Create combined variable for color mapping - to adapt
pcoa_df_clr$group <- paste0(pcoa_df_clr$control_or_infected, "_", pcoa_df_clr$light_or_dark)

pcoa_plot <- ggplot(pcoa_df_clr, aes(x = Axis1, y = Axis2, color = group, label = timepoint)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5, color = "black") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "control_light"   = "#7fbfff",  
      "control_dark"    = "#003f8c",  
      "infected_light"  = "#ff7f7f", 
      "infected_dark"   = "red"   
    ),
    name = "Condition"
  ) +
  labs(
    x = paste0("PCoA 1 (", round(pcoa_clr$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa_clr$values$Relative_eig[2] * 100, 1), "%)")
  )

ggsave("EDA/PCoA.png", pcoa_plot, width = 6, height = 5)