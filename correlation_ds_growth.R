# Compute and save correlation between defense systems expression and
# - bacterial growth (ribosomal genes expression)
# - phage growth (structural phage genes expression)


# Determine states to split (ex. infected and control, WT and mutant, stationary and exponential)
get_conditions <- function(metadata) {
  conds <- unique(metadata$State)
  conds <- conds[!is.na(conds)]
  
  if (length(conds) == 0) {
    return(NA)
  } else {
    return(conds)
  }
}


# correlation_df function
# 1. Compute the median clr of the proteins of the given defense system
# 2. Create a dataframe with 2 columns: the median clr of the defense system & the median clr of the marker
correlation_df <- function(proteins, markers_clr, condition = NA, use_all = FALSE) {

  # select clr for specified defense system protein
  if (use_all) {
    protein_clr <- clr_bacteria   
  } else {
    protein_clr <- clr[proteins, ]
  }
  
  # select column for given condition (always takes time 0)
  if (!is.na(condition)) {
    protein_clr <- protein_clr[,(str_detect(colnames(protein_clr), condition) | str_detect(colnames(protein_clr), paste0("^0(_|$)")))]
  }
  
  # compute median clr
  protein_clr <- protein_clr %>%
    summarise(across(everything(), median, na.rm = TRUE)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "run",
      values_to = "median_ds"
    )
  
  # Alignement of protein_clr and markers_clr
  corr_df <- protein_clr %>%
    left_join(markers_clr, by = "run")
  
  return(corr_df)
}

# compute_correlation function
# call correlation_df and return cor.test for specific ds and marker (+ the number of samples)
compute_correlation <- function(proteins, markers_clr, condition = NA, use_all = FALSE) {

  cor_df <- correlation_df(proteins, markers_clr, condition, use_all)
  n_samples <- nrow(cor_df)
  corr <- cor.test(cor_df$median_ds, cor_df$median_marker)
  return(list(corr = corr, n_samples = n_samples))
}
  
# compute_one_result
# call compute_correlation and returns a dataframe with all the info to store
compute_one_result <- function(ds_type, markers_clr, markers, condition = NA, use_all = FALSE) {
  
  if (use_all) {
    res <- compute_correlation(
      proteins = NULL,
      markers_clr = markers_clr,
      condition = condition,
      use_all = TRUE
    )
  } else {
    proteins <- ds %>% filter(type == ds_type) %>% pull(protein_in_syst)
    
    res <- compute_correlation(
      proteins = proteins,
      markers_clr = markers_clr,
      condition = condition,
      use_all = FALSE
    )
  }
  
  corr <- res$corr
  n_samples <- res$n_samples
  
  data.frame(
    Project = project,
    Species = species,
    Strain = strain,
    Phage = phage,
    Defense_system = ifelse(use_all, "ALL_GENES", ds_type),
    Condition = condition,
    marker_type = markers,
    corr = as.numeric(corr$estimate),
    pvalue = corr$p.value,
    n_samples = n_samples,
    stringsAsFactors = FALSE
  )
}



write_results_to_sheet <- function(df, sheet_name = "ds_correlations") {

  for (i in seq_len(nrow(df))) {
    sheet <- read_sheet(google_sheet_url, sheet = sheet_name)
    sheet_df <- as.data.frame(sheet)
    
    row <- df[i, ]

    row_id <- sheet_df %>%
      mutate(row_number = row_number()) %>%
      filter(
        Project == row$Project,
        Species == row$Species,
        Strain == row$Strain,
        Phage == row$Phage,
        Defense_system == row$Defense_system,
        (Condition == row$Condition |
           (is.na(Condition) & is.na(row$Condition)))
      ) %>%
      pull(row_number)

    # Column mapping
    col_corr  <- ifelse(row$marker_type == "ribosome", "G", "H")
    col_pval  <- ifelse(row$marker_type == "ribosome", "I", "J")
    col_padj  <- ifelse(row$marker_type == "ribosome", "K", "L")
    col_n     <- "S"
    
    # row already exists
    if (length(row_id) >= 1) {
      
      range_write(google_sheet_url, data.frame(row$corr),  range = paste0(col_corr, row_id + 1),  col_names = FALSE, sheet = sheet_name)
      range_write(google_sheet_url, data.frame(row$pvalue),range = paste0(col_pval, row_id + 1),  col_names = FALSE, sheet = sheet_name)
      range_write(google_sheet_url, data.frame(row$padj),  range = paste0(col_padj, row_id + 1),  col_names = FALSE, sheet = sheet_name)
      range_write(google_sheet_url, data.frame(row$n_samples), range = paste0(col_n, row_id + 1), col_names = FALSE, sheet = sheet_name)
      
    } else {
      
      new_row <- data.frame(
        Project = row$Project,
        Species = row$Species,
        Strain = row$Strain,
        Phage = row$Phage,
        Defense_system = row$Defense_system,
        Condition = row$Condition,
        Corr_ribosome = ifelse(row$marker_type == "ribosome", row$corr, NA),
        Corr_phage    = ifelse(row$marker_type == "phage", row$corr, NA),
        pvalue_ribosome = ifelse(row$marker_type == "ribosome", row$pvalue, NA),
        pvalue_phage    = ifelse(row$marker_type == "phage", row$pvalue, NA),
        padj_ribosome = ifelse(row$marker_type == "ribosome", row$padj, NA),
        padj_phage    = ifelse(row$marker_type == "phage", row$padj, NA),
        n_samples = row$n_samples
      )
      
      sheet_append(google_sheet_url, new_row, sheet = sheet_name)
    }
  }
}

conditions <- get_conditions(SRRs)

for (condition in conditions) {
  
  message("Running correlations for condition: ", condition)
  
  all_ds <- unique(ds$type)
  
  all_results <- purrr::map_dfr(all_ds, function(ds_type) {
    bind_rows(
      compute_one_result(ds_type, median_growth_markers_clr,  "ribosome", condition),
      compute_one_result(ds_type, median_struct_phage_clr, "phage", condition)
    )
  })
  
  # Add ALL_GENES
  all_results <- bind_rows(
    all_results,
    compute_one_result("ALL_GENES", median_growth_markers_clr,  "ribosome", condition, use_all = TRUE),
    compute_one_result("ALL_GENES", median_struct_phage_clr, "phage", condition, use_all = TRUE)
  )
  
  # padj
  all_results <- all_results %>%
    group_by(marker_type) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    ungroup()
  
  # write to sheet
  write_results_to_sheet(all_results)
}



