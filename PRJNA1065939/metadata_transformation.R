SRRs <- SRRs %>%
  mutate(time_point = as.numeric(str_remove_all(time_point, " min")), 
         genotype = str_replace_all(genotype, "rpoC G17D", "mutant")) %>%
  arrange(time_point)

SRRs <- SRRs %>%
  group_by(time_point, genotype) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(time_point, "_", genotype, "_", rep), State = genotype)