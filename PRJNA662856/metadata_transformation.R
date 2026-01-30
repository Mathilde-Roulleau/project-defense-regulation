SRRs <- SRRs %>%
  separate(Library.Name, into=c("Strain", "State", "rep"), sep = "_", remove=FALSE)%>%
  mutate(State = if_else(State == "uc", "control", "infected")) %>%
  mutate(condition = paste0(State, "_", rep), State = NA)