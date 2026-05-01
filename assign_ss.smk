import os
import glob
import subprocess
from pathlib import Path

SS_PREDICTOR_BINARY = os.environ.get("SS_PREDICTOR_BINARY", "stride")
PROTEIN_DIR = os.environ.get("PROTEIN_DIR", "protein_structures")

MOLECULE_IDS = [os.path.basename(f).replace(".pdb.gz", "") for f in glob.glob(f"{PROTEIN_DIR}/*.pdb.gz")]

rule all:
    input:
        expand("secondary_structure_data/{pdb}.ss.out", pdb=MOLECULE_IDS)

rule unzip_pdb:
    input:
        pdb = f"{PROTEIN_DIR}/{{pdb}}.pdb.gz"
    output:
        unzipped = temp("decompressed_files/{pdb}.pdb")
    shell:
        "zcat {input.pdb} > {output.unzipped}"

rule run_stride:
    input:
        unzipped = "decompressed_files/{pdb}.pdb"
    output:
        ss_out = "secondary_structure_data/{pdb}.ss.out"
    run:
        result = subprocess.run(
            [SS_PREDICTOR_BINARY, input.unzipped],
            capture_output=True, text=True
        )
        Path(output.ss_out).parent.mkdir(parents=True, exist_ok=True)
        Path(output.ss_out).write_text(result.stdout)
