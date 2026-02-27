### **Defense System Regulation in Bacteria Infected by Phages**

In this project, we are interested in understanding how expression of defense system genes is regulated in bacteria infected by phages. 

44 projects containing transcriptomic data from bacteria infected by a phage were selected from the NCBI database. After quality control steps 29 of them were used for the analysis. All of them are described in this [google sheet](https://docs.google.com/spreadsheets/d/1DwDgeSOugfCSfyAb6feXC-KmUsX-ejwrghj_iXVWOtg/edit?gid=0#gid=0). 

- **counts_merged.csv** – count data of all the projects obtained from alignment, counting and annotation

- Loading of reference genomes, defense systems (**Defense Finder**) and phage genes annotations (**Genomad**): 

    - Loading of **NCBI reference genomes** for all projects (stored in /project_id/ncbi_reference_genome folder) – cds_from_genomic.fna, genomic.gff, protein.faa
    
    - **order_proteins.py** – python script ordering proteins from protein.faa in the same order as genomic.gff and generate ordered_proteins.faa
    
    - **defense_finder.sh** – shell script (to be run with Linux) iterating on a list of projects to run Defense finder on ordered_proteins.faa (output defense_systems.tsv stored in /project_id/defense_finder folder)
    
    - **add_gene_id.py** – python script adding gene id to defense_systems.tsv (output defense_systems.csv stored in /project_id/defense_finder folder)
    
    - **get_phage_genes.sh** – shell script (to be run with Linux) iterating on a list of reference phage genomes (GCF_list.csv) to run Genomad (output reference_genome_genes.tsv stored in /project_id/genomad folder)

- **General_pipeline_iteration.ps1** – power shell script iterating on the given list of projects running the following scripts:

    - **General_pipeline.R** – call the following R scripts:
    
        - **load_packages.R** – load needed packages
        
        - **load_files.R** – load files of the project
        
        - **metadata_transformation.R** – specific to each project
        
        - **count_matrices.R** - create count matrixes for all genes + bacterial genes+ phage genes
        
        - **qc_steps.R** - quality control steps (percentage of detected genes and number of reads per run)
        
        - **CLR.R** – compute CLR from count matrices
        
        - **CLR_control.R** – plot CLR ranked 
        
        - **PCoA.R** – plot PCoA from CLR values (specific to each project)
        
        - **correlation_ds_growth.R** – compute correlation between DS genes and growth/infection (save them in the google sheet)
        
        - **partial_correlation_ds_growth.R** – compute partial correlation between DS genes and growth/infection (save them in the google sheet)
        
        - **correlations_all_genes.R** – compute correlation (simple and partial) between all the genes and growth/infection (save them in Correlation_all_genes_growth.csv)

- **DS_correlations_analysis** : 

    - ds_correlation_analysis.R – R script to plot results on defense system correlations
    
    - result plots

- **DS_partial_correlations_analysis**:

    - partial_correlation_analysis.R – R script to plot results on defense system partial correlations
    
    - result plots

- **All_genes_correlations_analysis**: 

    - All_genes_correlation_analysis.R - R script to plot results on all genes (simple and partial) correlations
    
    - result plots
