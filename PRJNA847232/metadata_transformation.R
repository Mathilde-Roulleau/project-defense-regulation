SRRs <- SRRs %>%
  group_by(plasmid) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(State = if_else(plasmid == "pNW129", "pNW129", "MotB"), condition = paste0(State, "_", rep))
