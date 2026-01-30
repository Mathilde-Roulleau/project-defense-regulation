SRRs <- SRRs %>%
  mutate(Time = as.numeric(if_else(
    str_detect(time_point_post_infection, "(Negative control)"), "0", time_point_post_infection))
  ) 

SRRs <- SRRs %>%
  group_by(Time) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(Time, "_", rep), State = NA)