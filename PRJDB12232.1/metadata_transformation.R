SRRs <- SRRs %>%
  mutate(Time = as.integer(sub("16_(\\d+)M(\\d+).*", "\\1", Sample_name)), 
         rep = sub("16_(\\d+)M(\\d+).*", "\\2", Sample_name))
SRRs <- SRRs %>% mutate(condition = paste0(time, "_", rep), State = NA)
