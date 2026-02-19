SRRs <- SRRs %>%
  mutate(Time = str_remove_all(infection_time, " min = no phage"))%>%
  mutate(Time = str_remove_all(Time, " min"))%>%
  mutate(rep = str_remove_all(Replicate_number, "experiment "))%>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(condition = paste0(Time, "_", genotype, "_", rep), State = genotype)
