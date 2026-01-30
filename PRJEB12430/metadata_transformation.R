SRRs <- SRRs %>%
  mutate(State = if_else(str_detect(Sample_name, "C"), "control", "infected"), 
         Time = as.numeric(if_else(str_detect(Sample_name, "8"), "8", "16")))

SRRs <- SRRs %>%
  group_by(Time, State) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", State, "_", rep))