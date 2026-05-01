#!/bin/bash
CORES=16
export PROTEIN_DIR="${1:-protein_structures}"
export SS_PREDICTOR_BINARY="${2:-/usr/local/bin/stride}"
echo " PDB to STRIDE"
snakemake -s assign_ss.smk --cores $CORES --keep-going
echo "STRIDE to Contexts"
snakemake -s build_contexts.smk --cores $CORES --keep-going
echo "Contexts to Angle list"
snakemake -s compute_angles.smk --cores $CORES --keep-going
echo "Plot"
Rscript analysis_scripts/create_density_plot.R
echo "Done. Output: outputs/orientation_density.png"
