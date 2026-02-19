# _____Load packages and google sheet containing correlation coefficients_____
source("load_packages.R")
google_sheet_url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1293325500#gid=1293325500"
p_correlation_ds_df <- read_sheet(google_sheet_url, sheet = "ds_correlations") %>%
  filter(no_sequencing_issue == "x", n_samples > 5)


# _____Plot the density correlation violin plot_____
df_filt <- p_correlation_ds_df %>%
  filter(Condition != "control" | is.na(Condition)) %>%
  mutate(density_axis = "density_axis")


tooltip_text <- function(corr_col) {
  aes(
    text = paste0(
      "Species: ", Species,
      "<br>Phage: ", Phage,
      "<br>Defense system: ", Defense_system,
      "<br>Correlation: ", round({{ corr_col }}, 3)
    )
  )
}

plot_corr_violin <- function(df, corr_col, ylab, color_by = NULL, orientation = c("horizontal", "vertical")) {
  
  orientation <- match.arg(orientation)
  
  if (orientation == "horizontal") {
    p <- ggplot(df) +
      geom_violin(
        aes(y = density_axis, x = {{ corr_col }}),
        color = "grey70", fill = "grey90", alpha = .2
      ) +
      geom_jitter(
        aes(y = density_axis, x = {{ corr_col }}),
            color = "grey40",
        width = 0, height = 0.2, size = 0.5, na.rm = TRUE
      ) +
      geom_vline(xintercept = 0, color = "grey") +
      #xlim(-1, 1) +
      xlab(ylab) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        panel.grid.minor = element_blank()
      )
  } else {
    p <- ggplot(df) +
      geom_violin(
        aes(x = "Structural phage", y = {{ corr_col }}),
        color = "grey70", fill = "grey90", alpha = .2
      ) +
      geom_jitter(
        aes(x = "Structural phage", y = {{ corr_col }},
            colour = {{ color_by }}),
        width = 0.2, height = 0, size = 0.5, na.rm = TRUE
      ) +
      geom_hline(yintercept = 0, color = "grey") +
      ylim(-1, 1) +
      ylab(ylab) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        panel.grid.minor = element_blank()
      )
  }
  
  p
}

p_ribo <- plot_corr_violin(df_filt,
  corr_col = partial_corr_ribosome,
  ylab = "Correlation to ribosomal genes independantly of infection effets",
  orientation = "horizontal")

ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/partial_correlations_analysis/density_partial_corr_ribosome.pdf",
  p_ribo, height = 50, width = 90, unit = "mm")

p_phage <- plot_corr_violin(df_filt,
  corr_col = partial_corr_phage, 
  ylab = "Correlation to structural phage genes independantly of growth effets",
  orientation = "horizontal")

ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/partial_correlations_analysis/density_partial_corr_phage.pdf",
  p_phage, height = 50, width = 90, unit = "mm")

# Highlighting projects
p_ribo_proj <- plot_corr_violin(df_filt, 
  partial_corr_ribosome,
  "Correlation to ribosomal genes",
  color_by = Project,
  orientation = "vertical")
p_phage_proj <- plot_corr_violin(
  partial_corr_phage,
  "Correlation to structural phage genes",
  color_by = Project,
  orientation = "vertical")

ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_ribosome_Xproject.pdf",
  p_ribo_proj, height = 90, width = 140, unit = "mm")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_phage_Xproject.pdf",
  p_phage_proj, height = 90, width = 140, unit = "mm")

#ggplotly(p_ribo_proj, tooltip = "text")

# Highlighting defense systems
p_ribo_ds <- plot_corr_violin(df_filt,
  partial_corr_ribosome,
  "Correlation to ribosomal genes",
  color_by = Defense_system,
  orientation = "vertical"
  )
p_phage_ds <- plot_corr_violin(df_filt,
  partial_corr_phage,
  "Correlation to structural phage genes",
  color_by = Defense_system,
  orientation = "vertical"
)

ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_ribosome_Xds.pdf",
  p_ribosome_ds, height = 90, width = 200, unit = "mm")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_phage_Xds.pdf",
  p_phage_ds, height = 90, width = 200, unit = "mm")

#ggplotly(p_phage_ds, tooltip = "text")



# _____Range correlation rho into categories_____
categorize_rho <- function(rho) {

  cut(rho,
      breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c(
        "very strong\u2212", "strong\u2212", "moderate\u2212", "weak\u2212", "very weak",
        "weak+", "moderate+", "strong+", "very strong+"), 
      ordered_result = TRUE)
}

p_correlation_ds_df <- p_correlation_ds_df %>%
  mutate(partial_corr_phage_strength    = categorize_rho(partial_corr_phage),
    partial_corr_ribosome_strength = categorize_rho(partial_corr_ribosome))


# _____Heat map of correlation strengths_____
heatmap_correlation_strength <- function(ds_or_project = c("project", "defense_system"),
                                         ribosome_or_phage = c("ribosome", "phage"),
                                         infected_or_control = c("infected", "control")) {
  
  ds_or_project <- match.arg(ds_or_project)
  ribosome_or_phage <- match.arg(ribosome_or_phage)
  
  # Select correlation columns
  corr_value <- if (ribosome_or_phage == "phage") "partial_corr_phage" else "partial_corr_ribosome"
  corr_strength <- if (ribosome_or_phage == "phage") {
    "partial_corr_phage_strength"
  } else {
    "partial_corr_ribosome_strength"
  }
  
  if (infected_or_control == "infected") {
    df <- p_correlation_ds_df %>%
      filter(Condition != "control" | is.na(Condition))
  } else {
    df <- p_correlation_ds_df %>%
      filter(Condition == "control")
  }
  
  if (ds_or_project == "project") {
    
    # Order projects by number of defense systems
    project_order <- p_correlation_ds_df %>%
      distinct(Project, Defense_system) %>%  
      group_by(Project) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n)) %>%
      pull(Project)
    
    df <- df %>%
      mutate(Project = factor(Project, levels = project_order)) %>%
      group_by(Project) %>%
      arrange(.data[[corr_value]], .by_group = TRUE) %>%
      mutate(tile_order = row_number()) %>%
      ungroup()
    
    p <- ggplot(df,
                aes(y = Project,
                    x = tile_order,
                    fill = .data[[corr_strength]])) +
      geom_tile(width = 0.9, height = 0.9, color = "white") +
      geom_text(data = df %>% filter(Defense_system != "ALL_GENES"),
        aes(label = Defense_system), size = 2, angle = 45) +
      geom_label(data = df %>% filter(Defense_system == "ALL_GENES"),
        aes(label = Defense_system), size = 2, angle = 45, fill = NA, linewidth = 0.4, color = "black") +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 7),
        panel.grid = element_blank()
      )
    
  } else {
    
    # Order defense systems by number of projects
    ds_order <- df %>%
      distinct(Defense_system, Project) %>%  
      group_by(Defense_system) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n)) %>%
      pull(Defense_system)
    
    df <- df %>%
      mutate(Defense_system = factor(Defense_system, levels = ds_order)) %>%
      group_by(Defense_system) %>%
      arrange(.data[[corr_value]], .by_group = TRUE) %>%
      mutate(tile_order = row_number()) %>%
      ungroup()
    
    p <- ggplot(df,
                aes(y = Defense_system,
                    x = tile_order,
                    fill = .data[[corr_strength]])) +
      geom_tile(width = 0.9, height = 0.9, color = "white") +
      geom_text(aes(label = Project),
                size = 0.7, angle = 45) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 7),
        panel.grid = element_blank()
      )
  }
  
  colors <- c(
    "very strong−" = "#313695",
    "strong−"      = "#4575b4",
    "moderate−"    = "#74add1",
    "weak−"        = "#abd9e9",
    "very weak"    = "#ffffbf",
    "weak+"        = "#fdae61",
    "moderate+"    = "#f46d43",
    "strong+"      = "#d73027",
    "very strong+" = "#a50026"
  )
  
  p +  scale_fill_manual(values = colors, name = paste(
        ifelse(ribosome_or_phage == "phage", "Phage", "Ribosome"),
        "correlation strength")) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9, face = "bold"),
      legend.key.size = unit(0.5, "cm")
    )
}


ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_partial_corr_ribosome_Xproject.pdf", 
       heatmap_correlation_strength("project", "ribosome", "infected"), height = 200, width = 200, unit = "mm")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_partial_corr_phage_Xproject.pdf", 
       heatmap_correlation_strength("project", "phage", "infected"), height = 200, width = 200, unit = "mm")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_partial_corr_ribosome_Xds.pdf", 
       heatmap_correlation_strength("defense_system", "ribosome", "infected"), height = 200, width = 200, unit = "mm")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_partial_corr_ribosome_Xds_control.pdf", 
       heatmap_correlation_strength("defense_system", "ribosome", "control"), height = 200, width = 200, unit = "mm")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/heatmap_partial_corr_ribosome_Xproject_control.pdf", 
       heatmap_correlation_strength("project", "ribosome", "control"), height = 200, width = 200, unit = "mm")



# _____Plot phage correlation vs ribosome correlation______
p_ribosome_vs_phage_ds_plot <- ggplot(p_correlation_ds_df %>% 
  filter(Defense_system != "ALL_GENES", Condition != "control" | is.na(Condition))) +
  geom_density_2d_filled(
    aes(x = partial_corr_ribosome, y = partial_corr_phage, fill = after_stat(level)),
    alpha = 0.4,
    contour = TRUE
  ) +
  scale_fill_grey(start = 1, end = 0, guide = "none") +
  geom_point(aes(x = Corr_ribosome, y = Corr_phage), color = "grey30") +
  #geom_point(aes(x = partial_corr_ribosome, y = partial_corr_phage, text = paste0("Species: ", Species,"<br>Phage: ", Phage,"<br>Defense system: ", Defense_system,"<br>Correlation_phage: ", round(partial_corr_phage, 3), "<br>Correlation_ribosome: ", round(partial_corr_ribosome, 3)), colour = Species),width = 0.2, height = 0,size = 0.5) +
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

p_ribosome_vs_phage_ds_plot <- ggMarginal(
  p_ribosome_vs_phage_ds_plot,
  type = "density",
  margins = "both",
  size = 5, color = "grey55",
  linewidth = 1.5
)

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/partial_correlations_analysis/ribosome_vs_phage_ds.pdf", 
       p_ribosome_vs_phage_ds_plot, height = 200, width = 200, unit = "mm")  
#ggplotly(p_ribosome_vs_phage_ds_plot, tooltip = "text", height = 800, width = 1200)

# _____Plot phage correlation vs ribosome correlation for ALL genes_____
ribosome_vs_phage_allgenes_plot <- ggplot(p_correlation_ds_df %>% 
  filter(Defense_system == "ALL_GENES", Condition != "control" | is.na(Condition))) +
  geom_density_2d_filled(
    aes(x = partial_corr_ribosome, y = partial_corr_phage, fill = after_stat(level)),
    alpha = 0.4,
    contour = TRUE
  ) +
  scale_fill_grey(start = 1, end = 0, guide = "none") +
  geom_point(aes(x = partial_corr_ribosome, y = partial_corr_phage), color = "grey20") +
  #geom_jitter(aes(x = partial_corr_ribosome, y = partial_corr_phage, text = paste0("Species: ", Species,"<br>Phage: ", Phage,"<br>Defense system: ", Defense_system,"<br>Correlation: ", round(partial_corr_ribosome, 3)), colour = Project),width = 0.2, height = 0,size = 0.5) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  ylab("Correlation to structural phage genes (all genes)") +
  xlab("Correlation to ribosomal genes (all genes)") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_vs_phage_all_genes.pdf", 
       ribosome_vs_phage_allgenes_plot, height = 200, width = 200, unit = "mm")  
#ggplotly(phage_corr_density_plot, tooltip = "text", height = 800, width = 500)




# _____Correlation rho vs Significance fdr padj_____
df_long <- p_correlation_ds_df %>%
  filter(Condition != "control" | is.na(Condition)) %>%
  select(Project, Defense_system,
         partial_corr_ribosome, partial_padj_ribosome,
         partial_corr_phage, partial_padj_phage) %>%
  pivot_longer(
    cols = c(partial_corr_ribosome, partial_corr_phage),
    names_to = "Type",
    values_to = "rho"
  ) %>%
  mutate(
    padj = ifelse(Type == "partial_corr_ribosome", partial_padj_ribosome, partial_padj_phage),
    Type = recode(Type,
                  partial_corr_ribosome = "Ribosome",
                  partial_corr_phage    = "Phage"))

rho_vs_padj <- ggplot(df_long, aes(x = rho, y = -log(padj))) +

  stat_density_2d(aes(fill = Type), geom = "polygon", alpha = 0.25, contour = TRUE) +
  geom_point(aes(color = Type), size = 2, alpha = 0.6) +
  geom_point(data = df_long %>% filter(Defense_system == "ALL_GENES"),
      aes(color = Type), shape = 21, fill = "white", size = 2, stroke = 1.2) +
  
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  
  xlim(-1, 1) +
  xlab("Pearson correlation") +
  ylab("-log(FDR-adjusted p-value)") +
  
  scale_color_manual(values = c("Ribosome" = "#1F78B4",
                                "Phage"    = "#E31A1C")) +
  scale_fill_manual(values = c("Ribosome" = "#1F78B4",
                               "Phage"    = "#E31A1C")) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/rho_vs_padj.pdf", 
       rho_vs_padj, height = 200, width = 200, unit = "mm")  


# _____Compute correlation correlation relative to global transcription_____
all_genes_corr <- p_correlation_ds_df %>%
  filter(Defense_system == "ALL_GENES") %>%
  select(Project, Condition,
         partial_corr_ribosome_all_genes = partial_corr_ribosome,
         partial_corr_phage_all_genes    = partial_corr_phage)

p_correlation_ds_df_normalized <- p_correlation_ds_df %>%
  left_join(all_genes_corr, by = c("Project", "Condition")) %>%
  mutate(
    partial_corr_ribosome_normalized = partial_corr_ribosome - partial_corr_ribosome_all_genes,
    partial_corr_phage_normalized    = partial_corr_phage - partial_corr_phage_all_genes
  ) %>%
  filter(Condition != "control" | is.na(Condition)) %>%
  mutate(density_axis = "density_axis")


# _____Density plot for relative correlation to global transcription_____
density_corr_ribo_normalized <- plot_corr_violin(
  p_correlation_ds_df_normalized, 
  partial_corr_ribosome_normalized, 
  "Project-normalized correlation to ribosomal genes", 
  "horizontal")
density_partial_corr_phage_normalized <- plot_corr_violin(
  p_correlation_ds_df_normalized, 
  partial_corr_phage_normalized, 
  "Project-normalized correlation to structural phage genes", 
  "horizontal")

ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_ribosome_normalized.pdf",
  density_corr_ribo_normalized, height = 50, width = 90, unit = "mm")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/density_partial_corr_phage_normalized.pdf",
  density_partial_corr_phage_normalized, height = 50, width = 90, unit = "mm")


# _____Plot phage correlation vs ribosome correlation (relative to global transcription)______

ribosome_vs_phage_ds_normalized <- ggplot(p_correlation_ds_df_normalized) +
  geom_density_2d_filled(
    aes(x = partial_corr_ribosome_normalized, y = partial_corr_phage_normalized, fill = after_stat(level)),
    alpha = 0.4,
    contour = TRUE
  ) +
  scale_fill_grey(start = 1, end = 0, guide = "none") +
  geom_point(aes(x = partial_corr_ribosome, y = partial_corr_phage), color = "grey20") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  ylab("Correlation to structural phage genes (relative to global transcription)") +
  xlab("Correlation to ribosomal genes (relative to global transcription)") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/ribosome_vs_phage_ds_normalized.pdf", 
       ribosome_vs_phage_ds_normalized, height = 200, width = 200, unit = "mm") 

# _____Boxplot per defense system of relative correlation to global transcription_____
# Only projects with more than 5 defense systems
boxplot_corr_ribo_normalized <- ggplot(p_correlation_ds_df_normalized, 
                                       aes(x = fct_reorder(Defense_system, partial_corr_ribosome_normalized, .fun = median, na.rm = TRUE), 
                                           y = partial_corr_ribosome_normalized)) +
  geom_boxplot(outlier.alpha = 0.3, color = "grey30", fill = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Correlation to ribosomal genes (relative to global transcription)") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 8),
        axis.title.y = element_blank())

boxplot_partial_corr_phage_normalized <- ggplot(p_correlation_ds_df_normalized, 
                                        aes(x = fct_reorder(Defense_system, partial_corr_phage_normalized, .fun = median, na.rm = TRUE), 
                                            y = partial_corr_phage_normalized)) +
  geom_boxplot(outlier.alpha = 0.3, color = "grey30", fill = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Correlation to structural phage genes (relative to global transcription)") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 8),
        axis.title.y = element_blank())

ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/boxplot_corr_ribo_normalized.pdf", 
       boxplot_corr_ribo_normalized, height = 180, width = 120, unit = "mm")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/correlations_analysis/boxplot_partial_corr_phage_normalized.pdf", 
       boxplot_partial_corr_phage_normalized, height = 180, width = 120, unit = "mm")

