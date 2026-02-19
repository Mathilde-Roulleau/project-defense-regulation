
# _____Get conditions_____
get_conditions <- function(metadata) {
  conds <- unique(metadata$State)
  conds <- conds[!is.na(conds)]
  if(length(conds) == 0) return(NA)
  conds
}

# _____Function to compute correlation_____
compute_one_result_allgenes <- function(gene, markers_clr, markers, condition = NA) {
  gene_clr <- clr[gene, ] %>%
    pivot_longer(cols = everything(), names_to = "run", values_to = "clr_gene")
  if (!is.na(condition)) {
    gene_clr <- gene_clr[(str_detect(gene_clr$run, condition) | str_detect(gene_clr$run, "^0(_|$)")), ]
  }
  corr_df <- gene_clr %>% left_join(markers_clr, by = "run")
  corr_res <- cor.test(corr_df$clr_gene, corr_df$median_marker)
  n_samples <- nrow(corr_df)
  data.frame(
    Project = project, Species = species, Strain = strain, Phage = phage,
    Gene_id = gene, Condition = condition, marker_type = markers,
    corr = as.numeric(corr_res$estimate),
    pvalue = corr_res$p.value,
    n_samples = n_samples,
    stringsAsFactors = FALSE
  )
}

# _____Function to compute partial correlation_____
compute_one_result_partial_allgenes <- function(gene, markers_clr_y, markers_y, markers_clr_z, markers_z, condition = NA) {
  gene_clr <- clr[gene, ] %>%
    pivot_longer(cols = everything(), names_to = "run", values_to = "clr_gene")
  if(!is.na(condition)) {
    gene_clr <- gene_clr[(str_detect(gene_clr$run, condition) | str_detect(gene_clr$run, "^0(_|$)")), ]
  }
  pcorr_df <- gene_clr %>%
    left_join(markers_clr_y, by = "run") %>%
    left_join(markers_clr_z, by = "run")
  
  pcorr <- pcor.test(
    x = pcorr_df$clr_gene,
    y = pcorr_df$median_marker.x,
    z = pcorr_df$median_marker.y
  )
  
  data.frame(
    Project = project, Species = species, Strain = strain, Phage = phage,
    Gene_id = gene, Condition = condition, marker_type = markers_y,
    partial_corr = as.numeric(pcorr$estimate),
    partial_pval = pcorr$p.value,
    stringsAsFactors = FALSE
  )
}

# _____Append or update CSV_____
update_csv_results <- function(df_new, file_path) {
  id_cols <- c("Project","Species","Strain","Phage","Gene_id","Condition","marker_type")
  df_new[, (id_cols) := lapply(.SD, as.character), .SDcols = id_cols]
  
  if(file.exists(file_path)) {
    dt <- fread(file_path)
    dt[, (id_cols) := lapply(.SD, as.character), .SDcols = id_cols]
    
    # Update existing rows or append
    setkeyv(dt, id_cols)
    setkeyv(df_new, id_cols)
    
    # Update correlations
    dt[df_new, `:=`(
      corr = i.corr,
      pvalue = i.pvalue,
      padj = i.padj,
      partial_corr = i.partial_corr,
      partial_pval = i.partial_pval,
      n_samples = i.n_samples
    )]
    
    # Add new rows
    dt <- rbind(dt, df_new[!dt, on = id_cols], fill=TRUE)
    
  } else {
    dt <- df_new
  }
  
  fwrite(dt, file_path)
}

# _____MAIN LOOP_____
file_path <- "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/Correlation_all_genes_growth.csv"
conditions <- get_conditions(SRRs)
all_genes <- rownames(clr_bacteria)

for(cond in conditions) {
  message("Processing condition: ", cond)
  
  all_results <- purrr::map_dfr(all_genes, function(gene) {
    # Correlation
    corr_ribo <- compute_one_result_allgenes(gene, median_growth_markers_clr, "ribosome", cond)
    corr_phage <- compute_one_result_allgenes(gene, median_struct_phage_clr, "phage", cond)
    
    # Partial correlation
    p_ribo <- compute_one_result_partial_allgenes(gene, median_growth_markers_clr, "ribosome",
                                         median_struct_phage_clr, "phage", cond)
    p_phage <- compute_one_result_partial_allgenes(gene, median_struct_phage_clr, "phage",
                                          median_growth_markers_clr, "ribosome", cond)
    
    # Combine into one row per gene per marker
    bind_rows(
      corr_ribo %>%
        left_join(p_ribo %>% select(Gene_id, partial_corr, partial_pval), by="Gene_id"),
      corr_phage %>%
        left_join(p_phage %>% select(Gene_id, partial_corr, partial_pval), by="Gene_id")
    )
  })
  
  # Adjust p-values for multiple testing
  all_results <- all_results %>%
    group_by(marker_type) %>%
    mutate(padj = p.adjust(pvalue, method = "BH"),
           partial_padj = p.adjust(partial_pval, method="BH")) %>%
    ungroup()
  
  # Write/update CSV
  update_csv_results(as.data.table(all_results), file_path)
}
  
  