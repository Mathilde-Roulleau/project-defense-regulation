# _____Load packages and google sheet containing correlation coefficients_____
source("load_packages.R")
google_sheet_url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1293325500#gid=1293325500"
correlation_df <- read_sheet(google_sheet_url, sheet = "Growth_markers") %>%
  filter(no_sequencing_issue == "x", n_samples > 5)

# Number of projects available for analysis
n_projects <- length(unique(correlation_df$Project))

# Number of unique defense systems
n_ds <- length(unique(correlation_df$Defense_system))

# Number of unique strain and phage
n_strain <- length(unique(correlation_df$Strain))
n_phage <- length(unique(correlation_df$Phage))

# Stats of the correlations gathered
stats_ribosome_correlation <- data.frame(
  mean   = mean(correlation_df$Corr_ribosome, na.rm = TRUE),
  median = median(correlation_df$Corr_ribosome, na.rm = TRUE),
  sd     = sd(correlation_df$Corr_ribosome, na.rm = TRUE)
)

stats_phage_correlation <- data.frame(
  mean   = mean(correlation_df$Corr_phage, na.rm = TRUE),
  median = median(correlation_df$Corr_phage, na.rm = TRUE),
  sd     = sd(correlation_df$Corr_phage, na.rm = TRUE)
)


# _____Plot defense systems distribution_____
ds_dist_plot <- ggplot(
  correlation_df,
  aes(
    x = fct_infreq(Defense_system),  
    fill = after_stat(count)      
  )
) +
  geom_bar() +
  scale_fill_gradient(low = "red1", high = "red4") +
  theme_minimal() +
  ylab("Defense system occurance in all the samples") +
  xlab(NULL) +
  theme(
    panel.grid = element_blank(),                 
    axis.text.x = element_text(angle = 90, size = 3, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 4),
    axis.title.y = element_text(size = 3),
    legend.position = "none"
  )

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ds_distribution.pdf", height = 45, width = 80, unit = "mm")  


# _____Plot number of defense system per strain_____
nb_ds_per_strain_plot <- ggplot(correlation_df %>%
  group_by(Species, Strain) %>%
  summarise(n_unique_defense = n_distinct(Defense_system), .groups = 'drop')%>%
    mutate(label = paste(Species, Strain)),
  aes(
    x = label,
    y = n_unique_defense,
    fill = n_unique_defense
  )
) +
  geom_col() +
  scale_fill_gradient(low = "red1", high = "red4") +
  theme_minimal() +
  ylab("Number of defense systems") +
  xlab(NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, size = 3, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 4),
    axis.title.y = element_text(size = 3),
    legend.position = "none"
  )

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/nb_ds_per_strain_plot.pdf", height = 45, width = 50, unit = "mm")  


# _____Plot the density vs correlation violin plot_____
## Highlighting projects
# Ribosome 

ribosome_corr_density_plot <- ggplot(correlation_df) +
  geom_violin(
    aes(x = "Structural phage", y = Corr_ribosome),
    color = "grey70", fill = "grey90", alpha = .2
  ) +
  geom_jitter(
    aes(x = "Structural phage", y = Corr_ribosome,
        text = paste0("Species: ", Species,
                      "<br>Phage: ", Phage,
                      "<br>Defense system: ", Defense_system,
                      "<br>Correlation: ", round(Corr_ribosome, 3)), colour = Project),
    width = 0.2, height = 0,
    size = 0.5
  ) +
  geom_hline(yintercept = 0, color = "grey") +
  ylim(-1, 1) +
  ylab("Correlation to ribosomal genes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_correlation_density.pdf", height = 90, width = 140, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)

# Phage
phage_corr_density_plot <- ggplot(correlation_df) +
  geom_violin(
    aes(x = "Structural phage", y = Corr_phage),
    color = "grey70", fill = "grey90", alpha = .2
  ) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(
    aes(x = "Structural phage", y = Corr_phage,
        text = paste0("Species: ", Species,
                      "<br>Phage: ", Phage,
                      "<br>Defense system: ", Defense_system,
                      "<br>Correlation: ", round(Corr_phage, 3)), colour = Project),
    width = 0.2, height = 0,
    size = 0.5
  ) +
  ylim(-1, 1) +
  ylab("Correlation to structural phage genes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/phage_correlation_density.pdf", height = 90, width = 140, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)


## Highlighting defense systems
# Ribosome

ribosome_corr_density_plot <- ggplot(correlation_df) +
  geom_violin(
    aes(x = "Structural phage", y = Corr_ribosome),
    color = "grey70", fill = "grey90", alpha = .2
  ) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(
    aes(x = "Structural phage", y = Corr_ribosome,
        text = paste0("Species: ", Species,
                      "<br>Phage: ", Phage,
                      "<br>Defense system: ", Defense_system,
                      "<br>Correlation: ", round(Corr_ribosome, 3)), colour = Defense_system),
    width = 0.2, height = 0,
    size = 0.5
  ) +
  ylim(-1, 1) +
  ylab("Correlation to ribosomal genes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_correlation_density_bis.pdf", height = 90, width = 200, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)

# Phage
phage_corr_density_plot <- ggplot(correlation_df) +
  geom_violin(
    aes(x = "Structural phage", y = Corr_phage),
    color = "grey70", fill = "grey90", alpha = .2
  ) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(aes(x = "Structural phage", y = Corr_phage,
        text = paste0("Species: ", Species,
                      "<br>Phage: ", Phage,
                      "<br>Defense system: ", Defense_system,
                      "<br>Correlation: ", round(Corr_phage, 3)), colour = Defense_system),
    width = 0.2, height = 0,
    size = 0.5
  ) +
  ylim(-1, 1) +
  ylab("Correlation to structural phage genes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/phage_correlation_density_bis.pdf", height = 90, width = 200, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)


# _____Range correlation rho into categories_____
categorize_rho <- function(rho) {

  cut(rho,
      breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c(
        "very strong\u2212", "strong\u2212", "moderate\u2212", "weak\u2212", "very weak",
        "weak+", "moderate+", "strong+", "very strong+"))
}

correlation_df <- correlation_df %>%
  mutate(Corr_phage_strength    = categorize_rho(Corr_phage),
    Corr_ribosome_strength = categorize_rho(Corr_ribosome))


# _____Heat map of correlation strengths_____
# Order projects by number of defense systems
project_order <- correlation_df %>%
  distinct(Project, Defense_system) %>%  
  group_by(Project) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  pull(Project)


# Order defense systems within each project by correlation
correlation_df <- correlation_df %>%
  mutate(Project = factor(Project, levels = project_order)) %>%
    group_by(Project) %>%
  arrange(Corr_phage, .by_group = TRUE) %>%
  mutate(tile_order = row_number()) %>%
  ungroup()

colors <- brewer.pal(9, "RdYlBu")  

# Plot stacked tiles
heatmap_stack <- ggplot(correlation_df, aes(x = Project, y = tile_order, fill = Corr_phage_strength)) +
  geom_tile(width = 0.9, height = 0.9, color = "white") +  
  geom_text(aes(label = Defense_system), size = 2, color = "black", angle = 45) + 
  scale_fill_manual(values = colors, name = "Phage correlation strength") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.5, "cm")
  )

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_corr_strength.pdf", height = 200, width = 200, unit = "mm")  


# _____Plot phage correlation vs ribosome correlation______
ribosome_vs_phage_plot <- ggplot(correlation_df) +
  geom_density_2d_filled(
    aes(x = Corr_ribosome, y = Corr_phage, fill = after_stat(level)),
    alpha = 0.4,
    contour = TRUE
  ) +
  scale_fill_grey(start = 1, end = 0, guide = "none") +
  geom_point(aes(x = Corr_ribosome, y = Corr_phage), color = "grey") +
  #geom_jitter(aes(x = Corr_ribosome, y = Corr_phage, text = paste0("Species: ", Species,"<br>Phage: ", Phage,"<br>Defense system: ", Defense_system,"<br>Correlation: ", round(Corr_ribosome, 3)), colour = Project),width = 0.2, height = 0,size = 0.5) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  ylab("Correlation to structural phage genes") +
  xlab("Correlation to ribosomal genes") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_vs_phage.pdf", height = 200, width = 200, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)

# _____Plot phage correlation vs ribosome correlation_____
ribosome_vs_phage_plot <- ggplot(correlation_df) +
  geom_density_2d_filled(
    aes(x = Corr_ribosome, y = Corr_phage, fill = after_stat(level)),
    alpha = 0.4,
    contour = TRUE
  ) +
  scale_fill_grey(start = 1, end = 0, guide = "none") +
  geom_point(aes(x = Corr_ribosome, y = Corr_phage), color = "grey") +
  #geom_jitter(aes(x = Corr_ribosome, y = Corr_phage, text = paste0("Species: ", Species,"<br>Phage: ", Phage,"<br>Defense system: ", Defense_system,"<br>Correlation: ", round(Corr_ribosome, 3)), colour = Project),width = 0.2, height = 0,size = 0.5) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  ylab("Correlation to structural phage genes") +
  xlab("Correlation to ribosomal genes") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_vs_phage.pdf", height = 200, width = 200, unit = "mm")  
ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)
