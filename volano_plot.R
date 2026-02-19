# _____Plot Volcano plot control vs infected_____

volcano_plot <- function(cols, condition, time) {
  # construction of DESeq dataset
  coldata <- data.frame(
    row.names = cols,
    condition = factor(condition)
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
  res<- mutate(res, phage = row.names(res) %in% row.names(cnt_matrix_phage))
  
  # plot -Log10 FDR vs Log fold change
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = ds)) +
    geom_point(data = subset(res, ds == "FALSE"), color = "grey", alpha = 0.3) +
    geom_point(data = subset(res, ds == "TRUE"), color = "darkgreen") +
    geom_point(data = subset(res, phage == "TRUE"), color = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_label_repel(data = subset(res, ds == "TRUE" & padj < 0.05),
                     aes(label = label), 
                     size = 2, 
                     max.overlaps = Inf, 
                     show.legend = NA,
                     colour = "black"
    ) +
    xlim(-20, 20) +
    ylim(0, 50) +
    annotate("text", x=0, y=45, label = paste0(unique(condition)[1], " vs ", unique(condition)[2]," (time ", time, "h)")) +
    theme_minimal(base_size = 10)
  
  p
}

plots <- list()

for(time in unique(SRRs$Time)){
  #for(light in unique(SRRs$light_or_dark)){
    cols <- SRRs%>%
      filter(Time==time)%>%
      #filter(Time == time, light_or_dark == light) %>%
      pull(condition) 
    condition <- SRRs%>%
      filter(Time==time)%>%
      #filter(Time == time, light_or_dark == light) %>%
      pull(State)
    plots[[paste0(time, "_", light)]] <- volcano_plot(cols, condition, time)
 # }
}

ncol <- 3
n_plots <- length(plots)
nrow <- ceiling(n_plots / ncol)

combined_volcano_plots <- wrap_plots(plots, ncol = ncol)

width_per_col <- 3
height_per_row <- 4

ggsave("plots_volcano&CLR/c_vs_i.png", combined_volcano_plots, width  = ncol * width_per_col,
       height = nrow * height_per_row, limitsize = FALSE)
