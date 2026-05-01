import sys
import os
import gzip
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm

window_dir = sys.argv[1]
result_file = sys.argv[2]
target_residue = sys.argv[3]
protein_dir = sys.argv[4] if len(sys.argv) > 4 else "protein_structures"

RESIDUE_BULK_CLASS = {
    "K": "Large",  "R": "Bulky",        "D": "Intermediate", "E": "Large",
    "G": "Tiny",   "A": "Tiny",         "V": "Small",        "L": "Intermediate",
    "I": "Intermediate", "M": "Large",  "P": "Small",        "S": "Small",
    "T": "Small",  "N": "Intermediate", "Q": "Large",        "C": "Small",
    "F": "Bulky",  "Y": "Bulky",        "W": "Bulky",        "H": "Large"
}


def signed_angle(v1, v2, axis):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    axis = axis / np.linalg.norm(axis)
    return np.degrees(np.arctan2(np.dot(np.cross(v1, v2), axis), np.dot(v1, v2)))


def load_pdb(pdb_id):
    path = os.path.join(protein_dir, f"{pdb_id}.pdb.gz")
    if not os.path.exists(path):
        return None
    with gzip.open(path, "rt") as fh:
        return PDBParser(QUIET=True).get_structure(pdb_id, fh)


def iter_residue(structure, chain_id, resnum):
    try:
        return next(
            res for model in structure
            for res in model[chain_id]
            if res.id[1] == resnum
        )
    except (StopIteration, KeyError):
        return None


def ca_vector(res):
    return res["CA"].get_coord() if res and "CA" in res else None


def sc_centroid(res):
    backbone = {"N", "CA", "C", "O"}
    pts = [a.get_coord() for a in res if a.get_name() not in backbone]
    return np.mean(pts, axis=0) if pts else None


def measure_angles(contexts_dir):
    for fname in os.listdir(contexts_dir):
        if not fname.endswith(".tsv"):
            continue
        fpath = os.path.join(contexts_dir, fname)
        if os.path.getsize(fpath) == 0:
            continue
        try:
            df = pd.read_csv(fpath, sep="\t", header=None)
        except Exception:
            continue
        if df.empty:
            continue

        pdb_id = df.iat[0, 12]
        structure = load_pdb(pdb_id)
        if structure is None:
            continue

        for i in range(0, len(df) - 2, 3):
            try:
                pr, cr, nr = df.iloc[i], df.iloc[i + 1], df.iloc[i + 2]
            except IndexError:
                continue
            if cr.iat[0] != target_residue or cr.iat[11] != "HHH":
                continue

            chain = cr.iat[1]
            rp = iter_residue(structure, chain, int(pr.iat[2]))
            rc = iter_residue(structure, chain, int(cr.iat[2]))
            rn = iter_residue(structure, chain, int(nr.iat[2]))
            if not all([rp, rc, rn]):
                continue

            ca_p, ca_c, ca_n = ca_vector(rp), ca_vector(rc), ca_vector(rn)
            sc_p, sc_c = sc_centroid(rp), sc_centroid(rc)
            if any(v is None for v in [ca_p, ca_c, ca_n, sc_p, sc_c]):
                continue

            axis = ca_c - ca_p
            if np.linalg.norm(axis) == 0:
                continue
            try:
                ang = signed_angle(sc_p - ca_p, sc_c - ca_c, axis)
            except Exception:
                continue

            aa1 = pr.iat[9]
            sz = RESIDUE_BULK_CLASS.get(aa1, "Unknown")
            if sz != "Unknown":
                yield [pdb_id, aa1, sz, ang]


collected_measurements = list(tqdm(measure_angles(window_dir), desc="Processing tripeptide windows"))

os.makedirs(os.path.dirname(result_file), exist_ok=True)
pd.DataFrame(collected_measurements, columns=["pdb", "left_aa", "size_class", "angle"]).to_csv(
    result_file, sep="\t", index=False
)
print(f"Saved: {result_file}  |  Total angles: {len(collected_measurements)}")
