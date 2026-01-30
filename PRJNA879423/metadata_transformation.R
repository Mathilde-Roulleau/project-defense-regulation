SRRs <- SRRs %>%
  mutate(Time = str_replace_all(treatment, 'untreated', '0'))%>%
  mutate(Time = as.numeric(str_remove_all(Time, " min with phage")))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)