# Construction of a distance matrix for PCoA
dist_matrix_cts <- cts[, .(Run, Geneid, counts)] %>%
  group_by(Run, Geneid) %>%
  pivot_wider(
    names_from = Run,
    values_from = counts,
    values_fill = 0
  ) %>%
  as.data.frame()


rownames(dist_matrix_cts) <- dist_matrix_cts$Geneid
dist_matrix_cts <- dist_matrix_cts[ , -1]

dist_matrix_cts <- vegdist(t(dist_matrix_cts), method = "euclidean")

# Compute principal coordinate decomposition
pcoa <- pcoa(dist_matrix_cts)

# Plot PCoA2 vs PCoA2
pcoa_df <- data.frame(
  Run = rownames(pcoa$vectors),
  Axis1 = pcoa$vectors[, 1],
  Axis2 = pcoa$vectors[, 2]
)
pcoa_df <- left_join(pcoa_df, SRRs, by = "Run")

pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = time, label = time)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5, color = "black") +
  theme_minimal() +
  labs(
    x = paste0("PCoA 1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA 2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)"),
    color = "Condition")

ggsave("EDA/PCoA.png", pcoa_plot, width = 6, height = 5)