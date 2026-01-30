SRRs <- SRRs %>%
  mutate(Time = as.numeric(str_remove_all(Time, " min")))%>%
  arrange(Time)

SRRs <- SRRs %>%
  group_by(Time, cell_state) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", cell_state, "_", rep), State = cell_state)