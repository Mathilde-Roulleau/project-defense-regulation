## COMPUTES CLR FOR ALL GENES + RIBOSOMAL GENES + STRUCTURAL PHAGE GENES + BACTERIAL GENES + PHAGE GENES


# Add pseudocount to avoid counts = 0
pseudocount <- min(cnt_matrix[cnt_matrix > 0])/10 
# Transform the cts_matrix dataframe into a Deseq object
dds_for_vst <- cnt_matrix + pseudocount

dds_for_vst <- DESeqDataSetFromMatrix(round(dds_for_vst), as.data.frame(colnames(dds_for_vst)), design = ~ 1)
# Apply a variance stabilizing transformation (VST) to the count data
vst <- varianceStabilizingTransformation(dds_for_vst) %>%
  assay() %>%
  as.data.frame()


vst <- 2^(vst) # because vst function apply a log(2)


# Normalization by gene length
genes_length <- cts$Length[match(rownames(vst), cts$Geneid)]
vst <- vst/(genes_length * 10^(-3))


# Centered Log2 Ratio (CLR)
clr <- log2(vst / exp(colMeans(log(vst)))) 


# CLR for ribosomal genes
growth_markers_genes <- c("rplV", "rplD", "rplC", "rpsC", "rpsS", "rplB", "rplP", "rplE", "rplX", "rimM", "rplF", "rplN", "rpsH", "rpsE", "rpsQ", "rplR", "rpsD", "secY", "rplQ", "rplO", "rpsM", "rpoA", "rpsK", "rplA", "rplK", "rplM", "rplJ", "pnp", "mnmC", "rpsG", "tig", "nusG", "tsf", "rpsI", "truB", "purB", "gyrA", "rbfA", "spoU", "rpsB")

growth_markers_proteins <- bacterial_genome[bacterial_genome$gene %in% growth_markers_genes, ]$protein_id
growth_markers_proteins <- unlist(growth_markers_proteins[!is.na(growth_markers_proteins)]) 

growth_markers_clr <- clr[growth_markers_proteins, ]

# CLR structural phage genes
struct_phage_clr <- clr[unique(cts[cts$Gene_origin == "Phage" & cts$Start %in% structural_phage_genes$start, ]$Geneid), ]

# CLR bacterial genes
clr_bacteria <- clr[row.names(clr) %in% cts[cts$Gene_origin == "Bacteria", ]$Geneid,]

# CLR phage genes
clr_phage <- clr[row.names(clr) %in% cts[cts$Gene_origin == "Phage", ]$Geneid,]


# Cumputing the median clr of the growth markers 

median_growth_markers_clr <- growth_markers_clr %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "run",
    values_to = "median_marker"
  )

median_struct_phage_clr <- struct_phage_clr %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "run",
    values_to = "median_marker"
  )

