SRRs <- SRRs %>%
  group_by(timepoint, control_or_infected) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(timepoint, "_", control_or_infected, "_", rep)) %>%
  mutate(timepoint = as.numeric(timepoint), State = control_or_infected)