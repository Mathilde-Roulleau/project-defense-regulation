SRRs <- SRRs[SRRs$strain == "SH1000",]
SRRs <- SRRs %>%
  mutate(time_point = case_when(
    time == "before infection" ~ 0,
    TRUE ~ as.numeric(str_remove_all(time, " min"))), 
    State = case_when(
      time == "before infection" ~ "control",
      TRUE ~ "infected")
  ) %>%
  arrange(time_point)

SRRs <- SRRs %>%
  group_by(time, State) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(time_point, "_", State, "_", rep), State = NA)
