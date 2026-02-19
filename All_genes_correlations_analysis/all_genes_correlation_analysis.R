# _____Load packages and csv files containing correlation coefficients for all the genes_____
source("../load_packages.R")
correlation_df <- fread("Correlation_all_genes_growth.csv") %>%
  filter(n_samples > 5)

# Number of projects available for analysis
n_projects <- length(unique(correlation_df$Project))

# Number of unique strain and phage
n_strain <- length(unique(correlation_df$Strain))
n_phage <- length(unique(correlation_df$Phage))

correlation_df <- correlation_df %>%
  filter(Condition != "control" | is.na(Condition)) %>%
  mutate(density_axis = "density_axis")

# _____Plot the density correlation violin plot_____
plot_corr_density <- function(marker, corr, xlab) {
  
  corr_col <- if (corr == "simple") "corr" else "partial_corr"
  df <- correlation_df %>%
    filter(marker_type == marker)

    ggplot(df) +
    geom_violin(
      aes(y = density_axis, x = .data[[corr_col]]),
      color = "grey70", fill = "grey90", alpha = .2
    ) +
    geom_jitter(
      data = ~ dplyr::slice_sample(.x, prop = 0.005),  
      aes(y = density_axis, x = .data[[corr_col]]),
      color = "grey40",
      width = 0, height = 0.2, size = 0.2
    ) +
    geom_vline(xintercept = 0, color = "grey") +
    xlab(xlab) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      text = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
}



density_ribo <- plot_corr_density("ribosome", "simple", "Correlation of all genes to ribosomal genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/density_corr_ribosome.pdf",
  density_ribo, height = 50, width = 90, unit = "mm")

density_phage <- plot_corr_density("phage", "simple", "Correlation of all genes to structural phage genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/density_corr_phage.pdf",
  density_phage, height = 50, width = 90, unit = "mm")

density_partial_ribo <- plot_corr_density("ribosome", "partial", "Correlation of all genes to ribosomal genes independantly of infection effets")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/density_partial_corr_ribosome.pdf",
  density_partial_ribo, height = 50, width = 90, unit = "mm")

density_partial_phage <- plot_corr_density("phage", "partial", "Correlation of all genes to structural phage genes independantly of growth effets")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/density_partial_corr_phage.pdf",
  density_partial_phage, height = 50, width = 90, unit = "mm")


# ____Histogram distribution of correlation values____
compute_hist <- function(df, col, bins = 200) {
  x <- as.numeric(df[[col]])
  hist(x, breaks = bins, plot = FALSE)
}

plot_corr_hist <- function(marker = c("ribosome", "phage"),
                           corr = c("simple", "partial"),
                           xlab,
                           bins = 200) {
  
  marker_val <- match.arg(marker)
  corr_val   <- match.arg(corr)
  corr_col   <- if (corr_val == "simple") "corr" else "partial_corr"
  
  df <- correlation_df %>%
    filter(marker_type == marker_val) %>%
    mutate(corr_value = as.numeric(.data[[corr_col]]))
  
  h <- compute_hist(df, corr_col, bins)
  
  hist_df <- data.frame(
    x = h$mids,
    y = h$density
  )
  
    ggplot(hist_df, aes(x = x, y = y)) +
      geom_col(
        fill = "grey85",          
        width = diff(h$breaks)[1]
      ) +
      geom_density(
        data = df,
        aes(x = corr_value, y = ..density..),
        color = "grey30",       
        size = 1,               
        alpha = 0.7              
      ) +
    geom_vline(xintercept = 0, color = "grey") +
    ylim(0, 0.8) +
    xlab(xlab) +
    ylab("Density") +
    theme_minimal() +
    theme(
      text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

hist_ribo <- plot_corr_hist("ribosome", "simple", "Correlation of all genes to ribosomal genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/hist_corr_ribosome.pdf",
  hist_ribo, height = 50, width = 120, unit = "mm")

hist_phage <- plot_corr_hist("phage", "simple", "Correlation of all genes to structural phage genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/hist_corr_phage.pdf",
  hist_phage, height = 50, width = 120, unit = "mm")

hist_partial_ribo <- plot_corr_hist("ribosome", "partial", "Partial correlation of all genes to ribosomal genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/hist_partial_corr_ribosome.pdf",
  hist_partial_ribo, height = 50, width = 120, unit = "mm")

hist_partial_phage <- plot_corr_hist("phage", "partial", "Partial correlation of all genes to structural phage genes")
ggsave(
  "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/hist_partial_corr_phage.pdf",
  hist_partial_phage, height = 50, width = 120, unit = "mm")



# _____Plot phage correlation vs ribosome correlation______
ribosome_vs_phage <- function(correlation) {
  
  if(correlation == "simple") {
    df_wide <- correlation_df %>%
      filter(Condition != "control" | is.na(Condition)) %>%
      select(Project, Species, Strain, Phage, Gene_id, Condition, marker_type, corr) %>%
      pivot_wider(
        names_from = marker_type,
        values_from = corr
      )
  } else {
    df_wide <- correlation_df %>%
      filter(Condition != "control" | is.na(Condition)) %>%
      select(Project, Species, Strain, Phage, Gene_id, Condition, marker_type, partial_corr) %>%
      pivot_wider(
        names_from = marker_type,
        values_from = partial_corr
      )
  }
  df_wide <- df_wide %>%
    mutate(
      phage = as.numeric(phage),
      ribosome = as.numeric(ribosome)
    )
  
  p <- ggplot(df_wide, aes(x = ribosome, y = phage)) +
    geom_point(alpha = 0) +
    geom_density_2d_filled(
      aes(fill = after_stat(level)),
      alpha = 0.4,
      contour = TRUE
    ) +
    scale_fill_grey(start = 1, end = 0, guide = "none") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_hline(yintercept = 0, color = "grey") +
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1)) +
    ylab("Correlation of all genes to structural phage genes") +
    xlab("Correlation of all genes to ribosomal genes") +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  p <- ggMarginal(
    p,
    type = "density",
    margins = "both",
    size = 5,
    color = "grey55",
    linewidth = 1.5
  )
  
  p
}


plot_ribosome_vs_phage <- ribosome_vs_phage("simple")
plot_ribosome_vs_phage_partial <- ribosome_vs_phage("partial")
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/ribosome_vs_phage.pdf", 
       plot_ribosome_vs_phage, height = 200, width = 200, unit = "mm")  
ggsave("C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/All_genes_correlations_analysis/ribosome_vs_phage_partial.pdf", 
       plot_ribosome_vs_phage_partial, height = 200, width = 200, unit = "mm") 
#ggplotly(plot_ribosome_vs_phage, tooltip = "text", height = 800, width = 500)