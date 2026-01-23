cnt_matrix <- cts %>% 
  left_join(SRRs, by = "Run") %>%   
  select(Geneid, condition, counts) %>%   
  pivot_wider(names_from = condition,               
              values_from = counts,               
              values_fill = 0) %>%   
  as.data.frame() 

rownames(cnt_matrix) <- cnt_matrix$Geneid 
cnt_matrix <- cnt_matrix[ , -1] 
cnt_matrix <- cnt_matrix[, order(names(cnt_matrix))]


# count matrix of bacterial genes
cnt_matrix_bacteria <- cts[cts$Gene_origin == "Bacteria", ] %>% 
  # filter bacteria genes   
  left_join(SRRs, by = "Run") %>%   
  select(Geneid, condition, counts) %>%   
  pivot_wider(names_from = condition,               
              values_from = counts,               
              values_fill = list(counts = 0)) %>%   
  as.data.frame()  
rownames(cnt_matrix_bacteria) <- cnt_matrix_bacteria$Geneid 
cnt_matrix_bacteria <- cnt_matrix_bacteria[ , -1] 
cnt_matrix_bacteria <- cnt_matrix_bacteria[, order(names(cnt_matrix_bacteria))]



# count matrix of phage genes
cnt_matrix_phage <- cts[cts$Gene_origin == "Phage", ] %>% 
  # filter phage genes   
  left_join(SRRs, by = "Run") %>%   
  select(Geneid, condition, counts) %>%   
  pivot_wider(names_from = condition,               
              values_from = counts,               
              values_fill = list(counts = 0)) %>%   
  as.data.frame()  
rownames(cnt_matrix_phage) <- cnt_matrix_phage$Geneid 
cnt_matrix_phage <- cnt_matrix_phage[ , -1] 
cnt_matrix_phage <- cnt_matrix_phage[, order(names(cnt_matrix_phage))]