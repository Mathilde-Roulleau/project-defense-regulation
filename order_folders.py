import os
import shutil

# Dossier racine contenant tous tes projets (modifie si besoin)
BASE_DIR = r"C:\Users\mathi\OneDrive - epfl.ch\project-defense-regulation"


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

for project in os.listdir(BASE_DIR):
    project_path = os.path.join(BASE_DIR, project)

    if not os.path.isdir(project_path):
        continue

    print(f"\n📁 Processing project: {project}")

    # Create target directories
    dirs = {
        "ncbi_reference_genome": [],
        "EDA": [],
        "genomad": [],
        "defense_finder": [],
        "plots_ds_vs_growth": [],
        "plots_volcano&CLR": []
    }

    for d in dirs:
        ensure_dir(os.path.join(project_path, d))

    # Scan files in project folder
    for file in os.listdir(project_path):
        file_path = os.path.join(project_path, file)

        if not os.path.isfile(file_path):
            continue

        fname = file.lower()

        # Rules
        if any(x in fname for x in ["gff", "faa", "fna"]):
            shutil.move(file_path, os.path.join(project_path, "ncbi_reference_genome", file))

        elif file in ["reads_per_run.png", "detected_genes.png", "PCoA.png"]:
            shutil.move(file_path, os.path.join(project_path, "EDA", file))

        elif fname.endswith("_genes.tsv"):
            shutil.move(file_path, os.path.join(project_path, "genomad", file))

        elif "defense_systems" in fname:
            shutil.move(file_path, os.path.join(project_path, "defense_finder", file))

        elif "ds_vs" in fname:
            shutil.move(file_path, os.path.join(project_path, "plots_ds_vs_growth", file))

        elif fname.endswith(".png"):
            shutil.move(file_path, os.path.join(project_path, "plots_volcano&CLR", file))

    print("✅ Done")
