path <- getwd()  
project <- str_extract(path, "PRJ[^\\\\/]*")
google_sheet_url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=1293325500#gid=1293325500"

# load info on project
projects_list <- read_sheet(google_sheet_url, sheet = "projects_list")
species <- projects_list %>%
  filter(Project_ID == project) %>%
  pull(Species)

strain <- projects_list %>%
  filter(Project_ID == project) %>%
  pull(Strain) %>%
  unlist()

phage <- projects_list %>%
  filter(Project_ID == project) %>%
  pull(Phage)



# load metadata
SRRs <- read.csv("metadata.csv", header = TRUE)

# load gene counts for each run
cts <- fread("../counts_merged.tsv")[run_accession %in% SRRs$Run]
names(cts)[names(cts) == "run_accession"] <- "Run"

# select runs and classify gene origin
bacterial_genome <- read.gff("ncbi_reference_genome/genomic.gff")
bacterial_genome <- bacterial_genome %>%
  separate_rows(attributes, sep = ";") %>%        
  separate(attributes, into = c("key", "value"), sep = "=", fill = "right") %>%  
  pivot_wider(names_from = key, values_from = value)

cts <- cts %>%
  mutate(Gene_origin = ifelse(
    sapply(Chr, function(x){
      any(strsplit(x, ";")[[1]] %in% bacterial_genome$seqid)
    }),
    "Bacteria",
    "Phage"
  ))


cts$Geneid <- sub("^cds-", "", cts$Geneid)
cts$Geneid <- sub("^rna-", "", cts$Geneid)


# load defense systems of the bacteria obtained from defense finder
ds <- read.csv("defense_finder/defense_systems.csv", header = TRUE)[, c("Geneid", "protein_in_syst", "type")]


# recover phage accession (gcf) from google sheet
phage_accession <- read_sheet(google_sheet_url, sheet = "GCF_list")
phage_accession <- phage_accession[phage_accession$Project_ID == project, ]$phage_accession_for_analysis

# phage genome
phage_genome <- read_tsv(paste0("genomad/", phage_accession, "_genes.tsv"))
structural_phage_genes <- phage_genome %>%
  filter(replace_na(str_detect(annotation_description, "tail|capsid|head"), FALSE))
