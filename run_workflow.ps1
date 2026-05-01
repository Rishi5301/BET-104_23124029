# Windows PowerShell Pipeline for Tripeptide Analysis
# This script runs all four stages of the analysis

$CORES = 16
$PROTEIN_DIR = if ($args.Count -gt 0) { $args[0] } else { "protein_structures" }
$SS_PREDICTOR_BINARY = if ($args.Count -gt 1) { $args[1] } else { "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" }

# Set environment variables
$env:PROTEIN_DIR = $PROTEIN_DIR
$env:SS_PREDICTOR_BINARY = $SS_PREDICTOR_BINARY

Write-Host " PDB to STRIDE" -ForegroundColor Cyan
snakemake -s assign_ss.smk --cores $CORES --keep-going

Write-Host "STRIDE to Contexts" -ForegroundColor Cyan
snakemake -s build_contexts.smk --cores $CORES --keep-going

Write-Host "Contexts to Angle list" -ForegroundColor Cyan
snakemake -s compute_angles.smk --cores $CORES --keep-going

Write-Host "Plot" -ForegroundColor Cyan
$RscriptPath = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
if (Test-Path $RscriptPath) {
    & $RscriptPath "analysis_scripts/create_density_plot.R"
} else {
    Write-Host "Error: Rscript not found at $RscriptPath" -ForegroundColor Red
    exit 1
}

Write-Host "Done. Output: outputs/orientation_density.png" -ForegroundColor Green
