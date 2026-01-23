# Plot CLR for each defense systems

plot_CLR <- function(protein, ds_type) {
  protein_clr <- clr[protein, ]
  
  protein_clr <- protein_clr %>%
    t() %>%
    as.data.frame() %>%
    rename(value = protein) %>%
    rownames_to_column("time_rep") %>%
    separate(time_rep, into = c("time", "rep"), sep = "_") 
  
  # compute mean and se
  median_clr <- protein_clr %>%
    group_by(time) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE), 
      .groups = "keep"
    ) %>%
    mutate(time = factor(time, levels = c("Pre-Lysis", "Infection"))) %>%
    arrange(time)
  
  # Plot
  ggplot(median_clr, aes(x = time, y = median, group = 1)) +
    geom_errorbar(aes(ymin = median-sd, ymax = median+sd)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_minimal() +
    ylim(-3.5, 3.5) +
    labs(
      title = paste("Defense system = ", ds_type),
      y = "CLR"
    )
}

plots <- list()

for (i in rownames(ds)) {
  plots[[i]] <- plot_CLR(ds[i, ]$protein_in_syst, ds[i, ]$type)
}

wrap_plots(plots[1:6], ncol = 3)


###################

# DESEQ + VOLCANO PLOT FUNCTION (control vs infected)


# columns of cnt_matrix to select
cols <- c("Pre-Lysis_1", "Pre-Lysis_2", "Infection_1", "Infection_2")

# construction of DESeq dataset
coldata <- data.frame(
  row.names = cols,
  condition = factor(c("control", "control", "infected", "infected"))
)

countData <- round(cnt_matrix[,cols])

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ condition)

# differential expression analysis
dds <- DESeq(dds)
res <- as.data.frame(results(dds, alpha = 0.1, pAdjustMethod = "fdr"))

# label for defense genes 
labels <- ds$type[match(rownames(res), ds$protein_in_syst)]
labels[is.na(labels)] <- ""

# add column labels (character and logical)
res <- mutate(res, label=labels)
res <- mutate(res, ds = factor(res$label != ""))

# plot -Log10 FDR vs Log fold change
plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = ds)) +
  geom_point(data = subset(res, ds == "FALSE"), color = "grey", alpha = 0.3) +
  geom_point(data = subset(res, ds == "TRUE"), color = "darkgreen") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_label_repel(data = subset(res, ds == "TRUE" & padj < 0.05),
                   aes(label = label), 
                   size = 3, 
                   max.overlaps = Inf, 
                   show.legend = NA,
                   colour = "black"
  ) +
  xlim(-5, 5) +
  ylim(0, 15) +
  annotate("text", x=0, y=13, label = "C vs I") +
  theme_bw()


plot

###########

#FUNCTION DESEQ + VOLCANO PLOT (time 0 vs time t)

deseq_volcano_0_t <- function(t, state) {
  # columns of cnt_matrix to select
  if (state == "control") {
    cols <- c(
      paste0("0_", state, "_1"),
      paste0("0_", state, "_2"),
      paste0(t, "_", state, "_1"),
      paste0(t, "_", state, "_2"))
  } else if (t == 240) {
    cols <- c(
      paste0("0_", state, "_1"),
      paste0("0_", state, "_2"),
      paste0("0_", state, "_3"),
      paste0(t, "_", state, "_1"),
      paste0(t, "_", state, "_2"))
  } else {
    cols <- c(
      paste0("0_", state, "_1"),
      paste0("0_", state, "_2"),
      paste0("0_", state, "_3"),
      paste0(t, "_", state, "_1"),
      paste0(t, "_", state, "_2"),
      paste0(t, "_", state, "_3"))
  }
  
  
  # construction of DESeq dataset
  if (state == "control") {
    coldata <- data.frame(
      row.names = cols,
      condition = factor(c("0", "0", "t", "t"))
    )
  } else if (t == 240) {
    coldata <- data.frame(
      row.names = cols,
      condition = factor(c("0", "0", "0", "t", "t"))
    )
  } else {
    coldata <- data.frame(
      row.names = cols,
      condition = factor(c("0", "0", "0", "t", "t", "t"))
    )
  }
  
  countData <- round(cnt_matrix[,cols])
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = coldata,
                                design = ~ condition)
  
  # differential expression analysis
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, alpha = 0.3))
  
  # label for defense genes 
  labels <- rownames(res) %>%
    as.data.frame() %>%
    `colnames<-`("protein_in_syst") %>%
    left_join(ds, by = "protein_in_syst") %>%
    mutate(label = ifelse(is.na(type), "", type)) %>%
    pull(label)
  
  # add column labels (character and logical)
  res <- mutate(res, label=labels)
  res <- mutate(res, ds = factor(res$label != ""))
  
  # plot -Log10 FDR vs Log fold change
  plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = ds)) +
    geom_point(data = subset(res, ds == "FALSE"), color = "grey", alpha = 0.3) +
    geom_point(data = subset(res, ds == "TRUE"), color = "darkgreen") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_label_repel(data = subset(res, ds == "TRUE"& padj < 0.05),
                     aes(label = label), 
                     size = 3, 
                     max.overlaps = Inf, 
                     show.legend = NA,
                     colour = "black"
    ) +
    xlim(-12, 12) +
    ylim(0, 40) +
    annotate("text", x=-4, y=40, label = paste0(state, " - time 0 VS ", t)) +
    theme_bw()
  
  
  return(plot)
  
}

plot_c_240 <- deseq_volcano_0_t(240, "control")
plot_i_240 <- deseq_volcano_0_t(240, "infected")

plot_c_240 + plot_i_240 + plot_layout(nrow = 1) 



##############

#PLOT DEFENSE SYSTEMS CLR vs GROWTH MARKERS CLR

plot_correlation <- function(protein, ds_type, markers_clr, markers, condition_filter=NULL) {
  protein_clr <- clr[protein, , drop = FALSE]
  
  protein_clr <- protein_clr %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("time_rep") %>%
    separate(time_rep, into = c("time", "rep"), sep = "-") %>%
    mutate(time = as.numeric(time)) %>%
    pivot_longer(
      cols = -c(time, rep),   
      names_to = "protein",
      values_to = "value"
    ) %>%
    group_by(time, rep) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # linear regression
  model <- lm(median ~ markers_clr, data = protein_clr)
  slope <- coef(model)[2]            # coefficient
  intercept <- coef(model)[1]        # intercept
  slope_label <- paste0("Slope = ", round(slope, 3))
  
  # add slope coefficient to google sheet
  url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1994194302#gid=1994194302"
  sheet_name <- "Growth_markers"
  sheet <- read_sheet(url, sheet = sheet_name)
  sheet_df <- as.data.frame(sheet)
  
  # Check if list already exists
  row_id <- sheet_df %>%
    mutate(row_number = row_number()) %>%
    filter(Project == project,
           Defense_system == ds_type) %>%
    pull(row_number)
  
  # Determine which column to fill
  col_id <- ifelse(markers == "ribosome genes",
                   "D",
                   "E")
  
  if (length(row_id) == 1) {
    # if row already exists
    
    range <- paste0(col_id, row_id + 1) 
    range_write(
      ss = url,
      data = data.frame(slope),
      col_names = FALSE, 
      range = range,
      sheet = sheet_name
    )
    
  } else {
    # row do not exists yet
    
    new_row <- data.frame(
      Project = project,
      Defense_system = ds_type,
      Condition = NA, 
      Slope_ribosome = ifelse(markers == "ribosome genes", slope, NA),
      Slope_phage    = ifelse(markers == "structural phage genes", slope, NA)
    )
    
    sheet_append(
      ss = url,
      data = new_row,
      sheet = sheet_name
    )
  }
  
  
  # color palette for time
  n <- length(unique(protein_clr$time))  
  palette_time <- brewer.pal(n, "Spectral")
  
  # Plot
  plot <- ggplot(protein_clr, aes(x = markers_clr, 
                                  y = median, 
                                  color = as.factor(time))) +
    geom_point(size = 5) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    annotate("text",
             x = min(markers_clr),
             y = max(protein_clr$median),
             label = slope_label,
             hjust = 0, vjust = 1,
             size = 5) + 
    scale_color_manual(values = palette_time) +
    labs(
      x = paste0("Median ", markers, " CLR"),
      y = paste0(ds_type, " CLR"),
      color = "Time"
    ) +
    theme_minimal()
  
  return(plot)
}


plots <- list()
for (ds_ in unique(ds$type)) {
  proteins <- ds %>% filter(type == ds_) %>% pull(protein_in_syst)
  plots[[ds_]] <- plot_correlation(proteins, ds_, median_growth_markers_clr, "ribosome genes")
}

wrap_plots(plots, ncol = 3)


plots <- list()
for (ds_ in unique(ds$type)) {
  proteins <- ds %>% filter(type == ds_) %>% pull(protein_in_syst)
  plots[[ds_]] <- plot_correlation(proteins, ds_, median_struct_phage_clr, "structural phage genes")
}

wrap_plots(plots, ncol = 3)


