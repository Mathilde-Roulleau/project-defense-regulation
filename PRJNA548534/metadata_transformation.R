SRRs <- SRRs %>%
  mutate(Time = case_when(
    str_detect(Time, "control") ~ 0,
    TRUE ~ as.numeric(str_extract(Time, "\\d+"))
  )
  )

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)