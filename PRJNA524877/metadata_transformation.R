SRRs <- SRRs %>%
  mutate(State = ifelse(str_detect(Sample.Name, "Preinfection"), "control", "infected"),
         time_point = str_extract(Sample.Name, "\\d+")
  ) 
SRRs <- SRRs %>%
  mutate(time_point = ifelse(is.na(time_point), 0, time_point), 
         condition = paste0(time_point, "_", State), State = NA)