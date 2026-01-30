SRRs <- SRRs %>%
  mutate(State = ifelse(str_detect(strain, "P4"), "infected", "control"),
         Time = as.numeric(str_remove_all(Time_point, " min post infection")))

# metadata not very clear on which sample is infected/control. As there is phage gene in the same proportions in all the samples, we consider all the samples to be infected
# exception for one sample at 3.5min. We remove it.
SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)

SRRs <- SRRs[!(SRRs$condition == "3.5_2"), ]

cts <- cts %>%
  filter(Run %in% SRRs$Run)