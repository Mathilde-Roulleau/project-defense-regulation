SRRs <- SRRs %>%
  mutate(Time = str_remove_all(infection_time, " minutes"))%>%
  mutate(Time = as.numeric(str_remove_all(Time, " min")))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)