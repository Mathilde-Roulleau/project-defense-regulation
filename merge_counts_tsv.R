source("load_packages.R")

sheet_url <- "https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg"

projects <- read_sheet(sheet_url, sheet = "Status")

projects_overlap <- projects %>%
  filter(counts_batch1 == "x", counts_batch2 == "x") %>%
  pull(Project_ID)


base_dir <- "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation"

get_srrs <- function(project_id) {
  metadata_path <- file.path(base_dir, project_id, "metadata.csv")
  
  if (!file.exists(metadata_path)) {
    warning(paste("metadata.csv manquant pour", project_id))
    return(character(0))
  }
  
  read_csv(metadata_path, show_col_types = FALSE) %>%
    pull(Run)
}

srrs_to_remove <- map(projects_overlap, get_srrs) %>%
  unlist() %>%
  unique()


counts_batch0 <- read_tsv("counts_batch0.tsv", col_names = TRUE)
counts_batch1 <- read_tsv("counts_batch1.tsv", col_names = TRUE)
counts_batch2 <- read_tsv("counts_batch2.tsv", col_names = TRUE)

counts_batch1_filtered <- counts_batch1 %>%
  filter(!run_accession %in% srrs_to_remove)

counts_batch2 <- counts_batch2 %>%
  mutate(
    Start = as.numeric(Start),
    End   = as.numeric(End)
  )

counts_merged <- bind_rows(list(
  counts_batch0,
  counts_batch1_filtered,
  counts_batch2)
)

write_tsv(counts_merged, "counts_merged.tsv")


