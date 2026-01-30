SRRs <- SRRs %>%
  mutate(Time = str_remove_all(Time, "h"))%>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(infection = str_replace_all(infection, "Alderaan", "infected"))%>%
  mutate(infection = str_replace_all(infection, "uninfected", "control")) 


SRRs <- SRRs %>%
  group_by(Time, infection) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", infection, "_", rep), State = infection)