SRRs <- SRRs %>%
  separate(Library.Name, into=c("Time", "rep"), sep = "-", remove=FALSE)%>%
  mutate(Time = as.numeric(str_remove_all(Time, "PQ"))) %>%
  mutate(condition = paste0(Time, "_", rep), State = NA) 