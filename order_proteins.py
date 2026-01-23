#!/usr/bin/env python

import re
import pathlib
import argparse
import Bio.SeqIO.FastaIO as FastaIO

folders = [
    "C:/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJEB31642/Bacillus cereus/ncbi_reference_genome"
]


def order_proteins(faa_file: pathlib.Path, gff_file: pathlib.Path, output_file: pathlib.Path):

    proteins = dict()

    with open(faa_file, "r") as handle:
        for (header, sequence) in FastaIO.SimpleFastaParser(handle):
            head, sep, tail = header.partition(" ")
            proteins[head] = sequence 

    header_order = list()

    with open(gff_file, "r") as handle:
        for line in handle:
            if not line.startswith("#"):
                line_split = line.rstrip().split("\t")
                feature = line_split[2]
                if feature == "CDS":
                    info = line_split[8]
                    if info.startswith("ID="):
                        id = info.split(";")[3]
                        if id.startswith("Name="):
                            header_order.append(id.split("=")[1])
                        else:
                            print(line)

    with open(output_file, "w") as target:
        for id in header_order:
            target.write(f">{id}\n")
            target.write(proteins[id] + "\n")

if __name__ == "__main__":
    for folder in folders:
        folder_path = pathlib.Path(folder)
        faa_file = folder_path / "protein.faa"
        gff_file = folder_path / "genomic.gff"
        output_file = folder_path / "ordered_protein.faa"

        # Check files exist
        if not faa_file.exists():
            print(f"Missing file: {faa_file}")
            continue
        if not gff_file.exists():
            print(f"Missing file: {gff_file}")
            continue

        print(f"Processing {folder_path.name}...")
        order_proteins(faa_file, gff_file, output_file)