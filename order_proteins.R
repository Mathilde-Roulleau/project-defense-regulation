# Order proteins in the AA sequences file (.faa) based on the nucleotides sequences file (.fna)
# Necessitates packages Biostrings, tidyr, dplyr, stringr, seqinr

order_proteins <- function(faa_path, fna_path) {
  
  # LOADING .fna FILE AND EXTRACTING PROTEIN ID
  fna <- read.fasta(fna_path, as.string = TRUE) # extract .fna file (nucleotide sequences)
  
  fna <- data.frame(
    sequence = sapply(fna, function(x) paste(x, collapse = "")), # DNA sequence
    name = sapply(fna, function(x) attr(x, "name")), # name (lcl|NZ_MCTE02000001.1_cds_WP_408103861.1_1)
    annot = sapply(fna, function(x) attr(x, "Annot")), # Annot (>lcl|NZ_MCTE02000001.1_cds_WP_408103861.1_1 [locus_tag=BCV12_RS00005] [protein=hypothetical protein] [frame=3] [protein_id=WP_408103861.1] [location=<1..425] [gbkey=CDS])
    stringsAsFactors = FALSE
  ) 
  
  fna <- fna %>% mutate( 
    protein_id = str_extract(fna$annot, "(?<=protein_id=)[^]]+") # add column with only protein_id (WP_408103861.1)
  ) 
  row.names(fna) <- NULL # reset indexes
  
  fna <- fna[!is.na(fna$protein_id), ] # remove rows with no protein_id
  
  # LOADING .faa FILE TO REORDER
  faa <- read.fasta("protein.faa")
  faa_ids <- sub(" .*", "", names(faa))  # create list with protein_id
  
  # REORDERING .faa FILE FROM .fna FILE
  fna <- fna[fna$protein_id %in% faa_ids,]
  order_idx <- match(fna$protein_id, faa_ids) # match the order of the proteins from .fna
  ordered_faa <- faa[order_idx] # apply the right order to the .faa
  
  # SAVING THE ORDERED .faa FILE
  write.fasta(sequences = ordered_faa, names = sub("^>", "", sapply(ordered_faa, function(x) attr(x, "Annot"))),file.out = "reordered.faa")
}