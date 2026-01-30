SRRs <- SRRs %>%
  mutate(condition = str_remove_all(Library.Name, "FACHB-524_"))%>%
  mutate(condition = str_remove_all(condition, "infected_"))%>%
  separate(condition, into=c("Time", "rep"), sep = "h", remove=FALSE)%>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)