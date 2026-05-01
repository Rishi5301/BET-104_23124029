configfile: "config.yaml"
import glob
import os
import subprocess

TARGET_RESIDUE = config["target_residue"]
PROTEIN_IDS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("secondary_structure_data/*.ss.out")]

rule all:
    input:
        "outputs/measurement_angles.tsv"

rule calculate_angles:
    input:
        expand("tripeptide_windows/context_for_{aa}_in_{pdb}.tsv", pdb=PROTEIN_IDS, aa=[TARGET_RESIDUE])
    output:
        "outputs/measurement_angles.tsv"
    params:
        aa=TARGET_RESIDUE,
        struct_dir=os.environ.get("PROTEIN_DIR", "protein_structures")
    run:
        subprocess.run(
            ["python", "analysis_scripts/compute_orientation_angles.py",
             "tripeptide_windows", output[0], params.aa, params.struct_dir],
            check=True
        )
