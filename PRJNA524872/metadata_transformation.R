SRRs <- SRRs %>%
  mutate(State = ifelse(str_detect(Sample.Name, "Control"), "control", "infected"),
         Time = str_extract(Sample.Name, "\\d+"), 
         condition = paste0(Time, "_", State)
  )