SRRs <- SRRs %>%
  mutate(Time = as.integer(str_remove_all(time_point, "min")))
SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)
