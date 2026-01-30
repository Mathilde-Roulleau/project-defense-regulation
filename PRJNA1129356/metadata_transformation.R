SRRs <- SRRs %>%
  mutate(condition = str_extract(infection_timepoint, "\\d+"), 
         condition = ifelse(is.na(condition), 0, condition),
         State = NA)