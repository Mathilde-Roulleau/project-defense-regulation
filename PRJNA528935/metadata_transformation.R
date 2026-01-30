SRRs <- SRRs %>%
  mutate(Time = str_remove_all(treatment, "min"))%>%
  mutate(Time = str_remove_all(Time, "after phiYY phage infection of "))%>%
  mutate(Time = as.numeric(str_replace_all(Time, "without phiYY phage infection", "0")))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)