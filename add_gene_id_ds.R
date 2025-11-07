# Add Gene ID column to the defense systems file from the annotation features file (.gff)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("2 arguments for this script <ds_file> <gff_file>")
}

ds_file <- args[1]
gff_file <- args[2]

library(rtracklayer)
library(dplyr)
library(tidyr)

add_gene_id_ds(ds_file, gff_file)


add_gene_id_ds <- function(ds_file, gff_file) {
  # load annotation features .gff
  gff <- import(gff_file) 
  gff <- as.data.frame(gff, row.names = NULL)
  gff_cds <- gff[gff[, "type"] == "CDS", c("Name", "locus_tag")] # only keep coding genes, only keep the protein and gene ID
  colnames(gff_cds)[colnames(gff_cds) == "Name"] <- "protein_in_syst"
  
  # load defense systems of the bacteria obtained from defense finder
  defense_systems <- read.csv(ds_file, header = TRUE, sep = "")
  # separate rows with several associated proteins
  defense_systems <- defense_systems %>% separate_rows("protein_in_syst", sep = ",")
  
  # join with the protein ID
  defense_systems <- defense_systems %>%
    left_join(gff_cds, by = "protein_in_syst")
  names(defense_systems)[names(defense_systems)=="locus_tag"] <- "Geneid"
  
  # save the defense systems file
  write.csv(defense_systems, file = "defense_systems.csv")
}