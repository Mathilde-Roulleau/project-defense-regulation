plot_ds_vs_growth <- function(protein, ds_type, markers_clr, markers, condition_filter = NA) {
  # clr for specified defense system protein
  protein_clr <- clr[protein, , drop = FALSE]

  # compute the median clr of the proteins of the given defense system
  
  if (!is.na(condition_filter)) {
    # column of type time_condition1_condition2_rep
      if (grepl("^[^_]+_[^_]+$", condition_filter)) {
        conditions_id <- setdiff(names(markers_clr), c("time", "rep", "median"))
        conditions <- strsplit(condition_filter, "_")[[1]]

        protein_clr <- protein_clr %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column("time_condition1_condition2_rep") %>%
          separate(time_condition1_condition2_rep, into = c("time", "condition1", "condition2","rep"), sep = "_") %>%
          pivot_longer(
            cols = -c(time, condition1, condition2, rep),   
            names_to = "protein",
            values_to = "value"
          ) %>%
          group_by(time, condition1, condition2, rep) %>%
          summarise(
            median = median(as.numeric(value), na.rm = TRUE),
            .groups = "drop"
          )
        
        protein_clr <- protein_clr %>% filter((condition1 == conditions[1] & 
           condition2 == conditions[2])| time == 0) %>%
            rename(!!conditions_id[1] := condition1,
              !!conditions_id[2] := condition2)
          
        # Alignement of protein_clr and markers_clr
        plot_df <- protein_clr %>%
          left_join(
            markers_clr %>%
              select(time, conditions_id[1], conditions_id[2], rep, median_marker = median),
            by = c("time", conditions_id[1], conditions_id[2], "rep")
          )

      } else {
        # columns of type time_condition_rep
        protein_clr <- protein_clr %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column("time_condition_rep") %>%
          separate(time_condition_rep, into = c("time", "condition","rep"), sep = "_") %>%
          pivot_longer(
            cols = -c(time, condition, rep),   
            names_to = "protein",
            values_to = "value"
          ) %>%
          group_by(time, condition, rep) %>%
          summarise(
            median = median(as.numeric(value), na.rm = TRUE),
            .groups = "drop"
          )
        protein_clr <- protein_clr %>%
          filter(condition == condition_filter | time == 0)
        markers_clr <- markers_clr %>%
          filter(condition == condition_filter | time == 0)
        
        
        # Alignement of protein_clr and markers_clr
        plot_df <- protein_clr %>%
          left_join(
            markers_clr %>%
              select(time, condition, rep, median_marker = median),
            by = c("time", "condition", "rep")
          )
      }
  
    
  } else {
    # columns of type time_rep
    protein_clr <- protein_clr %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("time_rep") %>%
      separate(time_rep, into = c("time", "rep"), sep = "_") %>%
      pivot_longer(
        cols = -c(time, rep),   
        names_to = "protein",
        values_to = "value"
      ) %>%
      group_by(time, rep) %>%
      summarise(
        median = median(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(condition = NA)
    
    # Alignement of protein_clr and markers_clr
    plot_df <- protein_clr %>%
      left_join(
        markers_clr %>%
          select(time, rep, median_marker = median),
        by = c("time", "rep")
      )
    
  }

  
  # linear regression
  model <- lm(median ~ median_marker, data = plot_df)
  slope <- coef(model)[2]            # coefficient
  intercept <- coef(model)[1]        # intercept
  slope_label <- paste0("Slope = ", round(slope, 3))
  
  # add slope coefficient to google sheet
  url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1994194302#gid=1994194302"
  sheet_name <- "Growth_markers"
  sheet <- read_sheet(url, sheet = sheet_name)
  sheet_df <- as.data.frame(sheet)
  sheet_df <- sheet_df %>%
    mutate(Condition = na_if(Condition, ""))
  
  # Check if list already exists
  row_id <- sheet_df %>%
    mutate(row_number = row_number()) %>%
    filter(
      Project == project,
      Defense_system == ds_type,
      if (is.na(condition_filter)) is.na(Condition)
      else Condition == condition_filter
    ) %>%
    pull(row_number)
  
  # Determine which column to fill in the google sheet
  col_id <- ifelse(markers == "ribosome genes",
                   "D",
                   "E")
  
  # Store value in google sheet
  if (length(row_id) >= 1) {
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
      Condition = condition_filter, 
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
  plot <- ggplot(plot_df, aes(
    x = median_marker,
    y = median,
    color = as.factor(time)
  )) +
    geom_point(size = 5) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    annotate(
      "text",
      x = min(plot_df$median_marker, na.rm = TRUE),
      y = max(plot_df$median, na.rm = TRUE),
      label = slope_label,
      hjust = 0, vjust = 1,
      size = 5
    ) +
    scale_color_manual(values = palette_time) +
    labs(
      x = paste0("Median ", markers, " CLR"),
      y = paste0(ds_type, " CLR"),
      title = condition_filter,
      color = "Time"
    ) +
    theme_minimal()
  
  
  return(plot)
}