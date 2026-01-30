SRRs <- SRRs %>%
  mutate(Time = str_remove_all(source_name, "bacterial cells with phage infection_"))%>%
  mutate(Time = str_remove_all(Time, " min"))
SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)