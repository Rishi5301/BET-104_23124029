import sys
import re
import os
import pandas as pd

annotation_file = sys.argv[1]
result_file = sys.argv[2]
target_residue = sys.argv[3]

pdb_id = os.path.basename(annotation_file).split(".")[0]

AMINO_ACID_NOTATION = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def extract_annotations(filepath):
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith("ASG"):
                continue
            fields = line.split()
            if len(fields) < 10:
                continue
            m = re.match(r"\d+", fields[3])
            if not m:
                continue
            try:
                yield {
                    "res3": fields[1], "chain": fields[2],
                    "resnum": int(m.group()), "idx": int(fields[4]),
                    "ss_code": fields[5], "ss_name": fields[6],
                    "phi": float(fields[7]), "psi": float(fields[8]),
                    "area": float(fields[9]),
                }
            except (ValueError, IndexError):
                continue


def generate_windows(records):
    for i in range(1, len(records) - 1):
        yield records[i - 1], records[i], records[i + 1]


all_records = list(extract_annotations(annotation_file))

window_rows = []
for prev, center, nxt in generate_windows(all_records):
    if center["res3"] != target_residue:
        continue
    helix_seq = "".join(AMINO_ACID_NOTATION.get(r["res3"], "X") for r in (prev, center, nxt))
    helix_ss = prev["ss_code"] + center["ss_code"] + nxt["ss_code"]
    positions = (f"{center['chain']}:{prev['resnum']},"
                 f"{center['chain']}:{center['resnum']},"
                 f"{center['chain']}:{nxt['resnum']}")
    for res in (prev, center, nxt):
        window_rows.append([
            res["res3"], res["chain"], res["resnum"], res["idx"],
            res["ss_code"], res["ss_name"], res["phi"], res["psi"], res["area"],
            AMINO_ACID_NOTATION.get(res["res3"], "X"), helix_seq, helix_ss, pdb_id, positions
        ])

pd.DataFrame(window_rows).to_csv(result_file, sep="\t", index=False, header=False)
