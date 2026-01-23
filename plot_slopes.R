df <- read_sheet(
  "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg",
  sheet = "Growth_markers"
)

df <- df %>%
  mutate(
    X_combination = paste(Project, Species, Strain, Phage, Condition, sep = " | "),
    X_label = paste0(Species, "_", Condition)
  )

n_defense <- length(unique(df$Defense_system))
colors <- colorRampPalette(
  brewer.pal(11, "Spectral")
)(n_defense)

ggplot(df, aes(x = X_combination,
               y = Slope_ribosome,
               color = Defense_system)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(
    x = df$X_label,
    y = "Slope [CLR(defense_systems) vs CLR(ribosomes)]",
    color = "Defense system", size = 5
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    legend.position = "right"
  )

ggsave(
  filename = "slope_ribosome.png",
  plot = last_plot(),
  width = 45,
  height = 600,
  units = "in",
  dpi = 300
)

