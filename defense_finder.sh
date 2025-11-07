#!/bin/bash

# Liste des dossiers contenant ordered_protein.faa
folders=(
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1025116"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1065939"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1068994"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA340008"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA641380"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJEB12430"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA524872"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA633474"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA528935"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA534259"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA610075"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1230733"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1242552"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA1254618"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA481888"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA882106"
"/mnt/c/Users/mathi/OneDrive - epfl.ch/project-defense-regulation/PRJNA847232"
)

# Dossier temporaire pour travailler
tmp_dir=~/tmp_defense_finder
mkdir -p "$tmp_dir"

# Dossier pour stocker tous les résultats finaux
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