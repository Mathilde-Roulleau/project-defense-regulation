SRRs <- SRRs %>%
  mutate(Time = str_extract(Submitter_Id, "(?<=_t)\\d+"),
         peptide = ifelse(str_detect(Submitter_Id, "with_pep"), TRUE, FALSE)) 
SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)