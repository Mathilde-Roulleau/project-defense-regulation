SRRs <- SRRs %>%
  mutate(Time = case_when(
    treatment_time == "untreated" ~ 0,
    TRUE ~ as.numeric(str_remove_all(treatment_time, " min with phage")))) %>%
  arrange(Time)

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)