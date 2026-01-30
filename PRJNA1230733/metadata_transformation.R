SRRs <- SRRs %>%
  mutate(condition = str_remove_all(Library.Name, "Control_"))%>%
  separate(condition, into=c("Time", "rep"), sep = "h_Rep", remove=FALSE)%>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(condition = str_remove_all(condition, "Rep")) %>%
  mutate(condition = str_remove_all(condition, "h"), State = NA)