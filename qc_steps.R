# Number of bacterial and phage gene counts per run
cts_per_run_b_p <- cts %>%
  left_join(SRRs, by = "Run") # join with SRR to get the condition to label

cts_per_run_b_p <- cts_per_run_b_p %>%
  group_by(condition, Gene_origin) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  as.data.frame()

colnames(cts_per_run_b_p)[colnames(cts_per_run_b_p) == "Library.Name"] <- "Run"

# Plots bacteria and phage gene counts per run 
reads_per_run <-ggplot(cts_per_run_b_p, aes(x = condition, y = counts, fill = Gene_origin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Phage" = "red", "Bacteria" = "grey33")) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Run",
    y = "# reads",
    fill = "Gene origin"
  ) +
  theme(
    axis.text.x = element_text(angle = 90) 
  )

ggsave("EDA/reads_per_run.png")

# Compute percentage of detected genes
detected_bacterial_genes <- apply(cnt_matrix_bacteria, 2, function(x) {
  as.numeric(sum(x != 0) / nrow(cnt_matrix_bacteria) * 100)
})
detected_phage_genes <- apply(cnt_matrix_phage, 2, function(x) {
  as.numeric(sum(x != 0) / nrow(cnt_matrix_phage) * 100)
})

df_plot <- data.frame(
  Sample = colnames(cnt_matrix_bacteria),
  Bacteria = detected_bacterial_genes,
  Phage = detected_phage_genes
) %>%
  pivot_longer(cols = c(Bacteria, Phage),
               names_to = "Origin",
               values_to = "Percent")

detected_genes <-ggplot(df_plot, aes(x = Sample, y = Percent, fill = Origin)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("grey30", "red")) +
  labs(
    y = "Percentage of detected genes",
    x = NULL,
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top"
  )

ggsave("EDA/detected_genes.png")
