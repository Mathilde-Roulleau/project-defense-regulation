#!/usr/bin/env python3

import pathlib
import pandas as pd


folders = [
    "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA836150",
]



def add_gene_id_ds(ds_file, gff_file, output_file):
    """Merge Defense Finder TSV with GFF to add gene IDs."""

    # Load GFF file 
    gff = pd.read_csv(gff_file, sep="\t", comment="#", header=None, dtype=str)
    gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    # Filter for coding genes (CDS)
    gff_cds = gff[gff["type"] == "CDS"].copy()
    gff_cds["protein_in_syst"] = gff_cds["attributes"].str.extract(r"Name=([^;]+)")
    gff_cds["locus_tag"] = gff_cds["attributes"].str.extract(r"locus_tag=([^;]+)")
    gff_cds = gff_cds[["protein_in_syst", "locus_tag"]].dropna()

    # Load defense systems TSV obtained from defense finder
    defense_systems = pd.read_csv(ds_file, sep=None, engine="python") 
    defense_systems = defense_systems.assign(
        protein_in_syst=defense_systems["protein_in_syst"].str.split(",")
    ).explode("protein_in_syst")
    defense_systems["protein_in_syst"] = defense_systems["protein_in_syst"].str.strip()

    # Merge 
    defense_systems = defense_systems.merge(gff_cds, on="protein_in_syst", how="left")
    defense_systems.rename(columns={"locus_tag": "Geneid"}, inplace=True)

    # Save result
    defense_systems.to_csv(output_file, index=False)


if __name__ == "__main__":
    for folder in folders:
        folder_path = pathlib.Path(folder)

        ds_file = folder_path / "defense_systems.tsv"
        gff_file = folder_path / "genomic.gff"
        output_file = folder_path / "defense_systems.csv"

        if not ds_file.exists():
            print(f"Missing file: {ds_file}")
            continue
        if not gff_file.exists():
            print(f"Missing file: {gff_file}")
            continue

        add_gene_id_ds(ds_file, gff_file, output_file)
    
