#!/bin/bash

folders=(
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1065926"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1076276"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA644242"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA323606"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA345214"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJEB35542"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA662856"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA524877"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJEB31642"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA548534"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA530399"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA262612"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJDB12232"
    "/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1129356"
)

# Temporary folder
tmp_dir=~/tmp_defense_finder
mkdir -p "$tmp_dir"

# Folder to store results
results_dir=~/defense_finder_results
mkdir -p "$results_dir"

for folder in "${folders[@]}"; do
    echo "Processing $folder..."
    
    # Extraire le nom du projet pour renommer le fichier .tsv
    project_name=$(basename "$folder")
    
    # Copier le fichier ordered_protein.faa dans le dossier temporaire
    cp "$folder/ordered_protein.faa" "$tmp_dir/"
    
    # Se placer dans le dossier temporaire
    cd "$tmp_dir" || exit
    
    # Lancer defense-finder
    defense-finder run ordered_protein.faa
    
    # Renommer et déplacer le fichier .tsv vers le dossier final
    mv ordered_protein_defense_finder_systems.tsv "$results_dir/${project_name}_defense_finder_systems.tsv"
    
    # Supprimer les fichiers temporaires
    rm ordered*
    
    echo "Finished $folder."
done

# Ouvrir VSCode sur tous les résultats
code "$results_dir"

echo "All done! Results are in $results_dir"