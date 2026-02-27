#!/bin/bash

#############################
# ALLOWED PROJECT IDS
#############################

ALLOWED_PROJECTS=(
    "PRJNA323606.8"
    "PRJDB12232.4"
    "PRJNA430510.2"
    "PRJNA430510.3"
    "PRJNA430510.4"
)

is_allowed_project () {
    local pid="$1"
    for allowed in "${ALLOWED_PROJECTS[@]}"; do
        if [[ "$pid" == "$allowed" ]]; then
            return 0
        fi
    done
    return 1
}

### CONFIG ###
GCF_SHEET="/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/GCF_list.csv"
OUTPUT_BASE="/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation"
GENOMAD_DB="/root/genomad_db"
THREADS=8

#############################
# DETECT COLUMNS IN CSV
#############################

echo "[INFO] Detecting column numbers..."

col_PROJECT_ID=$(head -1 "$GCF_SHEET" | tr ',' '\n' | nl -v 1 | grep -w "Project_ID" | awk '{print $1}')
col_ACCESSION=$(head -1 "$GCF_SHEET" | tr ',' '\n' | nl -v 1 | grep -w "phage_accession_for_analysis" | awk '{print $1}')

if [[ -z "$col_PROJECT_ID" || -z "$col_ACCESSION" ]]; then
    echo "[ERROR] Missing required columns in CSV."
    exit 1
fi

echo "[INFO] Project_ID column = $col_PROJECT_ID"
echo "[INFO] phage_accession_for_analysis column = $col_ACCESSION"

echo "[INFO] Extracting required columns..."
csvcut -c "$col_PROJECT_ID","$col_ACCESSION" "$GCF_SHEET" > filtered.csv

#############################
# PROCESS ROWS
#############################

tail -n +2 filtered.csv | while IFS=',' read -r PROJECT_ID ACCESSIONS; do

    PROJECT_ID=$(echo "$PROJECT_ID" | tr -d '[:space:]')

    if [[ -z "$PROJECT_ID" || -z "$ACCESSIONS" ]]; then
        echo "[WARN] Empty or incomplete row — skipping."
        continue
    fi

    if ! is_allowed_project "$PROJECT_ID"; then
        echo "[INFO] Project $PROJECT_ID not in allowed list — skipping."
        continue
    fi


    echo "============================================"
    echo "[INFO] Processing PROJECT = $PROJECT_ID"
    echo "[INFO] Accessions = $ACCESSIONS"
    echo "============================================"

    TARGET_DIR="$OUTPUT_BASE/$PROJECT_ID"
    mkdir -p "$TARGET_DIR"

    #############################
    # SPLIT ACCESSIONS LIST
    #############################

    IFS=',' read -ra ACC_LIST <<< "$ACCESSIONS"

    for ACC in "${ACC_LIST[@]}"; do
        ACC=$(echo "$ACC" | tr -d '[:space:]')  # Trim spaces

        if [[ -z "$ACC" ]]; then
            echo "[WARN] Skipping empty accession entry."
            continue
        fi

        echo "----"
        echo "[INFO] Processing accession $ACC"
        echo "----"

        WORKDIR=$(mktemp -d)
        cd "$WORKDIR"

        echo "[INFO] Downloading genome $ACC"
        datasets download genome accession "$ACC" --include genome

        unzip -o ncbi_dataset.zip >/dev/null

        FASTA=$(find ncbi_dataset/data -type f -name "*genomic.fna" | head -n 1)

        if [[ ! -f "$FASTA" ]]; then
            echo "[ERROR] Genomic FASTA not found for $ACC"
            cd /
            rm -rf "$WORKDIR"
            continue
        fi

        echo "[INFO] Running genomad annotate on $ACC"
        genomad annotate --cleanup --splits "$THREADS" "$FASTA" genomad_output "$GENOMAD_DB"

        TSV=$(find genomad_output -type f -name "*genomic_genes.tsv" | head -n 1)

        if [[ ! -f "$TSV" ]]; then
            echo "[ERROR] genes.tsv NOT found for $ACC"
            cd /
            rm -rf "$WORKDIR"
            continue
        fi

        OUTFILE="$TARGET_DIR/${ACC}_genes.tsv"

        echo "[INFO] Saving → $OUTFILE"
        cp "$TSV" "$OUTFILE"

        cd /
        rm -rf "$WORKDIR"

    done

done

echo "[INFO] Pipeline complete!"

