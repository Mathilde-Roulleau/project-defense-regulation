#!/usr/bin/env python

import re
import pathlib
import argparse
import Bio.SeqIO.FastaIO as FastaIO

def get_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--faa', help='faa file of the proteins downloaded from ncbi.', required=True, type=pathlib.Path)
    parser.add_argument('-g', '--gff', help='gff file from ncbi.', required=True, type=pathlib.Path)
    parser.add_argument('-o', '--output', help='path and name of the output file.', required=True, type=pathlib.Path)
    return parser.parse_args()

args = get_arguments()
proteins = dict()

with open(args.faa, "r") as handle:
    for (header, sequence) in FastaIO.SimpleFastaParser(handle):
        head, sep, tail = header.partition(" ")
        proteins[head] = sequence 

header_order = list()

with open(args.gff, "r") as handle:
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

with open(args.output, "w") as target:
    for id in header_order:
        target.write(f">{id}\n")
        target.write(proteins[id] + "\n")