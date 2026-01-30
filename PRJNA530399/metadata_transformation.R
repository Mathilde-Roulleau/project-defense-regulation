SRRs <- SRRs %>%
  mutate(Time = coalesce(as.numeric(str_extract(source_name, "\\d+")), 0))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA
         )