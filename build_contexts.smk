import glob
import subprocess

configfile: "config.yaml"

TARGET_RESIDUE = config["target_residue"]
PROTEIN_IDS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("secondary_structure_data/*.ss.out")]

rule all:
    input:
        expand("tripeptide_windows/context_for_{aa}_in_{pdb}.tsv", pdb=PROTEIN_IDS, aa=[TARGET_RESIDUE])

rule extract_context:
    input:
        ss = "secondary_structure_data/{pdb}.ss.out"
    output:
        tsv = "tripeptide_windows/context_for_{aa}_in_{pdb}.tsv"
    params:
        aa = "{aa}"
    run:
        subprocess.run(
            ["python", "analysis_scripts/parse_secondary_structure.py", input.ss, output.tsv, params.aa],
            check=True
        )
