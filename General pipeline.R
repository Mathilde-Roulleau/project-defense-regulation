# General pipeline
# setwd() in console to set the desired project folder

# Load packages and files
#source("../load_packages.R")
source("../load_files.R")
source("metadata_transformation.R")

# Create count matrixes for all genes + bacterial genes+ phage genes
source("../count_matrices.R")

# Quality control steps
source("../qc_steps.R")

# CLR
source("../CLR.R")
source("../CLR_control.R")

# PCoA from CLR values
source("PCoA.R")

# Compute correlation between genes expression (ds and all genes) and bacterial growth / phage growth
source("../correlation_ds_growth.R")
#source("../plotCLR.R")
source("../partial_correlation_ds_growth.R")
source("../correlations_all_genes.R")
