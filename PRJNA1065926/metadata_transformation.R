SRRs <- SRRs %>%
  mutate(Time = c(15, 15, 15, 15, 30, 30, 30, 30))

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)
