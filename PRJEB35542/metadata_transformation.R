SRRs <- SRRs %>%
  mutate(Time = `Experimental_Factor._time..exp.`, 
         State = if_else(`Experimental_Factor._infect..exp.` == "none", "control", "infected"))%>%
  group_by(Time, State) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", State, "_", rep))