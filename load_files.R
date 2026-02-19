path <- getwd()  
project <- str_extract(path, "PRJ[^\\\\/]*")
google_sheet_url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1293325500#gid=1293325500"

# load info on project
projects_list <- read_sheet(google_sheet_url, sheet = "projects_list") %>%
  filter(Project_ID == project)

species <- projects_list$Species
strain  <- unlist(projects_list$Strain)
phage   <- projects_list$Phage


# load metadata
SRRs <- read.csv("metadata.csv", header = TRUE)

# load gene counts for each run
cts <- fread("../counts_merged.csv")[run_accession %in% SRRs$Run]
names(cts)[names(cts) == "run_accession"] <- "Run"

# select runs and classify gene origin
bacterial_genome <- as.data.table(read.gff("ncbi_reference_genome/genomic.gff"))

bacterial_genome[, gene :=
                   fifelse(
                     grepl("gene=", attributes),
                     sub(".*gene=([^;]+).*", "\\1", attributes),
                     NA_character_) ]
bacterial_genome[, protein_id :=
                   fifelse(
                     grepl("protein_id=", attributes),
                     sub(".*protein_id=([^;]+).*", "\\1", attributes),
                     NA_character_) ]

bacterial_genome <- bacterial_genome[, .(seqid, gene, protein_id)]

cts[, Gene_origin :=
      ifelse(vapply(strsplit(Chr, ";"),
                    function(x) any(x %in% bacterial_genome$seqid),
                    logical(1)),
             "Bacteria", "Phage")]



cts$Geneid <- sub("^cds-", "", cts$Geneid)
cts$Geneid <- sub("^rna-", "", cts$Geneid)


# load defense systems of the bacteria obtained from defense finder
ds <- read.csv("defense_finder/defense_systems.csv", header = TRUE)[, c("Geneid", "protein_in_syst", "type")]


# recover phage accession (gcf) from google sheet
phage_accession <- read_sheet(google_sheet_url, sheet = "GCF_list")
phage_accession <- phage_accession[phage_accession$Project_ID == project, ]$phage_accession_for_analysis

# phage genome
phage_genome <- read_tsv(paste0("genomad/", phage_accession, "_genes.tsv"))
structure_phage <- regex("tail|capsid|head", ignore_case = TRUE)
structural_phage_genes <- phage_genome %>%
  filter(str_detect(annotation_description, structure_phage))

