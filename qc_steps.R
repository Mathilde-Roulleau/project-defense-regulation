# Number of bacterial and phage gene counts per run
cts_per_run_b_p <- cts %>%
  left_join(SRRs, by = "Run") # join with SRR to get the condition to label

cts_per_run_b_p <- cts_per_run_b_p %>%
  group_by(condition, Gene_origin) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  as.data.frame()

colnames(cts_per_run_b_p)[colnames(cts_per_run_b_p) == "Library.Name"] <- "Run"

# Plots bacteria and phage gene counts per run 
p <-ggplot(cts_per_run_b_p, aes(x = condition, y = counts, fill = Gene_origin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Phage" = "red", "Bacteria" = "darkblue")) +
  theme_minimal() +
  labs(
    x = "Run",
    y = "# reads",
    fill = "Gene origin"
  ) +
  theme(
    axis.text.x = element_text(angle = 90) 
  )

print(p)

# Compute percentage of detected genes
detected_bacterial_genes <- apply(cnt_matrix_bacteria, 2, function(x) {
  as.numeric(sum(x != 0) / nrow(cnt_matrix_bacteria) * 100)
})
detected_phage_genes <- apply(cnt_matrix_phage, 2, function(x) {
  as.numeric(sum(x != 0) / nrow(cnt_matrix_phage) * 100)
})

mat <- rbind(detected_bacterial_genes, detected_phage_genes)

barplot(mat,
        beside = TRUE,   
        col = c("darkblue", "red"),
        ylab = "percentage of detected genes",
        ylim = c(0, 100), 
        las = 2, 
        cex.names = 0.5)

legend("bottomright", legend = c("Bacteria", "Phage"),
       fill = c("darkblue", "red"))