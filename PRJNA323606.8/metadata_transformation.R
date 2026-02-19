SRRs <- SRRs %>%
  mutate(condition = if_else(str_detect(infection_timepoint, "no exogenous phage added"), "control", "infected"),
         State = condition)