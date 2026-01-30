# Construction of a distance matrix for PCoA
dist_matrix_clr <- vegdist(t(clr), method = "euclidean")

# Compute principal coordinate decomposition
pcoa_clr <- pcoa(dist_matrix_clr)

# Plot PCoA2 vs PCoA2
pcoa_df_clr <- data.frame(
  condition = rownames(pcoa_clr$vectors),
  Axis1 = pcoa_clr$vectors[, 1],
  Axis2 = pcoa_clr$vectors[, 2]
)

pcoa_df_clr <- left_join(pcoa_df_clr, SRRs, by = "condition")

pcoa_plot <- ggplot(pcoa_df_clr, aes(x = Axis1, y = Axis2, color = State, label = State)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5, color = "black") +
  theme_minimal() +
  scale_color_manual(
    values = c("infected" = "red",  
               "control" = "darkblue")   
  ) +
  labs(
    x = paste0("PCoA 1 (", round(pcoa_clr$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa_clr$values$Relative_eig[2] * 100, 1), "%)"),
    color = "Condition")

ggsave("EDA/PCoA.png", pcoa_plot, width = 6, height = 5)