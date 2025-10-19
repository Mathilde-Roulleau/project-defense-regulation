# Add Gene ID column to the defense systems file from the annotation features file (.gff)
# Necessitates packages rtracklayer, dplyr

add_gene_id_ds <- function(ds_file, gff_file) {
  # load annotation features .gff
  gff <- import("genomic.gff") 
  gff <- as.data.frame(gff, row.names = NULL)
  gff_cds <- gff[gff[, "type"] == "CDS", c("Name", "locus_tag")] # only keep coding genes, only keep the protein and gene ID
  colnames(gff_cds)[colnames(gff_cds) == "Name"] <- "protein_in_syst"
  
  # load annotation features .gff
  gff <- import("genomic.gff") 
  gff <- as.data.frame(gff, row.names = NULL)
  gff_cds <- gff[gff[, "type"] == "CDS", c("Name", "locus_tag")] # only keep coding genes, only keep the protein and gene ID
  colnames(gff_cds)[colnames(gff_cds) == "Name"] <- "protein_in_syst"
  
  # join with the protein ID
  defense_systems <- defense_systems %>%
    left_join(gff_cds, by = "protein_in_syst")
  names(defense_systems)[names(defense_systems)=="locus_tag"] <- "Geneid"
  
  # save the defense systems file
  write.csv(defense_systems, file = "defense_systems.csv")
}