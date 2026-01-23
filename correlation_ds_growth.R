# Compute and save correlation between defense systems expression and
# - bacterial growth (ribosomal genes expression)
# - phage growth (structural phage genes expression)

correlation_df <- function(proteins, markers_clr, condition = NA) {
  # 1. Compute the median clr of the proteins of the given defense system
  # 2. Create a dataframe with 2 columns: the median clr of the defense system & the median clr of the marker

  # select clr for specified defense system protein
  protein_clr <- clr[proteins,]

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

compute_correlation <- function(proteins, markers_clr, condition = NA) {
  cor_df <- correlation_df(proteins, markers_clr, condition)
  n_samples <- nrow(cor_df)
  corr <- cor.test(cor_df$median_ds, cor_df$median_marker)
  return(list(corr = corr, n_samples = n_samples))
}
  
save_correlation <- function(ds_type, markers_clr, markers, condition = NA) {
  
  proteins <- ds %>% filter(type == ds_type) %>% pull(protein_in_syst)
  res <- compute_correlation(proteins, markers_clr, condition)
  corr <- res$corr
  n_samples <- res$n_samples
  
  corr_estimate <- corr$estimate
  corr_pvalue <- corr$p.value

  sheet_name <- "Growth_markers"
  sheet <- read_sheet(google_sheet_url, sheet = sheet_name)
  sheet_df <- as.data.frame(sheet)

  # Check if list already exists
  row_id <- sheet_df %>%
    mutate(row_number = row_number()) %>%
    filter(
      Project == project,
      Species == species,
      Strain == strain, 
      Phage == phage,
      Defense_system == ds_type,
      (Condition == condition |
          (is.na(Condition) & is.na(condition)))) %>%
    pull(row_number)
  
  # Determine which column to fill in the google sheet
  col_id_esti <- ifelse(markers == "ribosome",
                   "G",
                   "H")
  col_id_pvalue <- ifelse(markers == "ribosome",
                   "I",
                   "J")
  col_id_n <- "K"

    # Store value in google sheet
  if (length(row_id) >= 1) {
    # if row already exists
    range <- paste0(col_id_esti, row_id + 1) 
    range_write(
      ss = google_sheet_url,
      data = data.frame(corr_estimate),
      col_names = FALSE, 
      range = range,
      sheet = sheet_name
    )
    range <- paste0(col_id_pvalue, row_id + 1)
    range_write(
      ss = google_sheet_url,
      data = data.frame(corr_pvalue),
      col_names = FALSE, 
      range = range,
      sheet = sheet_name
    )
    range <- paste0(col_id_n, row_id + 1)
    range_write(
      ss = google_sheet_url,
      data = data.frame(n_samples),
      col_names = FALSE, 
      range = range,
      sheet = sheet_name
    )
    
  } else {
    # row do not exists yet
    new_row <- data.frame(
      Project = project,
      Species = species,
      Strain = strain, 
      Phage = phage,
      Defense_system = ds_type,
      Condition = condition, 
      Corr_ribosome = ifelse(markers == "ribosome", corr_estimate, NA),
      Corr_phage    = ifelse(markers == "phage", corr_estimate, NA),
      pvalue_ribosome = ifelse(markers == "ribosome", corr_pvalue, NA),
      pvalue_phage    = ifelse(markers == "phage", corr_pvalue, NA), 
      n_samples = n_samples
    )
    
    sheet_append(
      ss = google_sheet_url,
      data = new_row,
      sheet = sheet_name
    )
  }
}

