SRRs <- SRRs %>%
  mutate(condition = str_remove_all(Library.Name, "WC1-T"))%>%
  mutate(condition = str_remove_all(condition, "_R1_R2"))%>%
  separate(condition, into=c("Time", "rep"), sep = "-", remove=FALSE)%>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(condition = str_remove_all(condition, "Rep")) %>%
  mutate(condition = str_replace_all(condition, "-", "_"), State = NA)