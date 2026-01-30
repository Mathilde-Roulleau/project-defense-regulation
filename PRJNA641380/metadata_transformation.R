SRRs <- SRRs %>%
  mutate(TEX = ifelse(str_detect(source_name, "TEX"),
                      "TEX+", "TEX-"), 
         infection_status = str_replace_all(infection_status, " with phiKZ", "")) %>%
  mutate(infection_status = str_replace_all(infection_status, "uninfected", "control"))

SRRs <- SRRs %>%
  group_by(infection_status) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(condition = paste0(infection_status, "_", rep), State = infection_status) 
