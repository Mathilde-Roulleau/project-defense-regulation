SRRs <- SRRs%>% 
  rename(Time = stage_of_infection) %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)