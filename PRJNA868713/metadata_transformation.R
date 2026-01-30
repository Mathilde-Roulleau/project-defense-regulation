SRRs <- SRRs %>%
  mutate(Time = str_remove_all(time, " min"))%>%
  mutate(Time = as.numeric(Time))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)