# Base directory
$baseDir = "C:\Users\mathi\OneDrive - epfl.ch\project-defense-regulation"

# Path to your general R script
$rScript = Join-Path $baseDir "General pipeline.R"

# List of projects
$projects = @(
  "PRJEB35542","PRJNA1065939","PRJNA1068994","PRJNA1230733","PRJNA1242552","PRJNA252740",
  "PRJNA315575","PRJNA340008","PRJNA481888","PRJNA524872","PRJNA524877","PRJNA528935",
  "PRJNA530399","PRJNA534259","PRJNA548534","PRJNA610075","PRJNA662856","PRJNA681237",
  "PRJNA720911","PRJNA787900.1","PRJNA787900.2","PRJNA868713","PRJNA879423"
)

# Iterate over projects
foreach ($proj in $projects) {
    Write-Host "Processing project: $proj"
    
    # Go to project directory
    $projDir = Join-Path $baseDir $proj
    if (!(Test-Path $projDir)) {
        Write-Warning "Directory $projDir not found, skipping..."
        continue
    }
    Set-Location $projDir
    
    # Run R script
    Rscript $rScript
    
    Write-Host "Finished project: $proj"
}
