# Sidechain Orientation Analysis in Alpha-Helical Contexts

A computational pipeline for analyzing how the size of a preceding residue influences the rotameric state of arginine sidechains in protein alpha-helices.

---

## Research Question

In protein alpha-helices, arginine (ARG) sidechains can rotate freely within their rotameric wells. However, the immediately preceding residue (i−1) may physically constrain this rotation through steric effects. This pipeline investigates: **Does the size of residue i−1 predict how residue i (ARG) orients its sidechain?**

To answer this, the pipeline:
1. Identifies all alpha-helical ARG residues across the PDB
2. Computes the signed angle between consecutive CA→sidechain-centroid vectors
3. Groups results by the bulk classification of the preceding residue
4. Generates a density plot showing whether angle distributions differ across bulk classes

---

## Prerequisites & Installation

### System Requirements
- **Python 3.10+** with packages:
  - `biopython` — protein structure parsing
  - `numpy` — numerical computations
  - `pandas` — data manipulation
  - `tqdm` — progress bars

- **R 4.0+** with packages:
  - `ggplot2` — data visualization
  - `readr` — TSV file reading

- **External Tools:**
  - [STRIDE](http://webclu.bio.wzw.tum.de/stride/) — secondary structure prediction
  - [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0 — workflow management

### Setup

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Install R dependencies (if not already installed):**
   ```r
   install.packages(c("ggplot2", "readr"))
   ```

3. **Ensure STRIDE is accessible:**
   - Linux/macOS: Place or link `stride` to `/usr/local/bin/stride`
   - Windows: Update the `SS_PREDICTOR_BINARY` parameter in pipeline calls

---

## Quick Start

### Running the Pipeline

**Linux/macOS (Bash):**
```bash
chmod +x run_workflow.sh
./run_workflow.sh [PROTEIN_DIR] [SS_PREDICTOR_BINARY]
```

**Windows (PowerShell):**
```powershell
.\run_workflow.ps1
```

### Default Parameters
- `PROTEIN_DIR`: `protein_structures/` (directory containing `.pdb.gz` files)
- `SS_PREDICTOR_BINARY`: `/usr/local/bin/stride` (path to STRIDE executable)
- Target residue: `ARG` (configurable in `config.yaml`)

### Output Files
- `outputs/measurement_angles.tsv` — angle measurements with metadata
- `outputs/orientation_density.png` — density plot visualization
- `outputs/protein_catalog.txt` — list of processed PDB identifiers

---

## Project Structure

```
.
├── README.md                      This file
├── requirements.txt               Python dependencies (frozen versions)
├── config.yaml                    Configuration (target residue)
│
├── run_workflow.sh                Pipeline entry point (Bash)
├── run_workflow.ps1               Pipeline entry point (PowerShell)
│
├── assign_ss.smk                  Stage 1: Secondary structure assignment
├── build_contexts.smk             Stage 2: Tripeptide window extraction
├── compute_angles.smk             Stage 3: Angle calculation
│
├── analysis_scripts/              Analysis modules
│   ├── parse_secondary_structure.py     STRIDE output parser
│   ├── compute_orientation_angles.py    3D coordinate analysis
│   └── create_density_plot.R            Visualization renderer
│
├── protein_structures/            Input: Compressed PDB files (*.pdb.gz)
├── decompressed_files/            Temporary: Uncompressed PDBs
├── secondary_structure_data/      Intermediate: STRIDE annotations
├── tripeptide_windows/            Intermediate: Context windows (TSV)
│
└── outputs/                       Final results
    ├── measurement_angles.tsv           All computed angles
    ├── orientation_density.png          Colored density plot
    └── protein_catalog.txt              PDB identifiers processed
```

---

## Pipeline Overview

The analysis is divided into four sequential stages:

### Stage 1: Secondary Structure Assignment
**Input:** PDB structures (.pdb.gz)  
**Tool:** STRIDE  
**Output:** Secondary structure annotations

For each PDB file:
1. Decompress to temporary location
2. Run STRIDE secondary structure prediction
3. Save annotations to `secondary_structure_data/{pdb}.ss.out`
4. Clean up temporary files

### Stage 2: Tripeptide Context Extraction
**Input:** STRIDE annotations  
**Script:** `parse_secondary_structure.py`  
**Output:** Tripeptide windows (TSV format)

For each PDB's STRIDE output:
1. Parse all residue records using generator-based streaming
2. Extract sliding tripeptide windows (prev, center, next)
3. Filter for windows where center residue matches target (ARG)
4. Save context windows to `tripeptide_windows/context_for_{aa}_in_{pdb}.tsv`

Each output TSV contains three rows per tripeptide with 14 fields:
- **Structural:** residue name (3-letter), chain ID, residue number, DSSP index
- **Secondary structure:** SS code, SS name, phi/psi angles, accessible surface area
- **Context:** single-letter code, helix sequence, helix SS pattern, PDB ID, positions

### Stage 3: Angle Computation
**Input:** Tripeptide windows + PDB coordinates  
**Script:** `compute_orientation_angles.py`  
**Output:** Angle measurements (TSV format)

For each tripeptide window:
1. Load 3D coordinates from original PDB using BioPython
2. Verify all three residues are helical (HHH pattern required)
3. Calculate sidechain centroid vectors
4. Compute signed dihedral angle relative to helix axis
5. Assign preceding residue to bulk class
6. Append to results table

**Angle Definition:**
```
Vectors:
  sc_prev = centroid(prev sidechain) − CA(prev)
  sc_cent = centroid(ARG sidechain) − CA(ARG)
  axis    = CA(ARG) − CA(prev)

Angle:
  θ = arctan2((sc_prev × sc_cent) · axis, sc_prev · sc_cent)
```
Angle range: −180° to +180°. Positive indicates counter-clockwise rotation when viewed from prev toward ARG along helix axis.

### Stage 4: Visualization
**Input:** Angle measurements (TSV)  
**Script:** `create_density_plot.R`  
**Output:** Density plot (PNG)

Generate publication-quality density plot:
1. Read angle measurements
2. Normalize angles to −180° to +180° range
3. Assign colors to bulk classes (blue → green → orange → red palette)
4. Plot kernel density estimates for each class
5. Save as PNG with 300 DPI

---

## Residue Classification

The preceding residue is classified into one of five bulk categories based on sidechain size:

| Class | Residues | Sidechain Characteristics |
|-------|----------|--------------------------|
| **Tiny** | G, A | Minimal sidechain (H or CH₃) |
| **Small** | V, P, S, T, C | Short, compact (<4 non-H atoms) |
| **Intermediate** | L, I, N, D | Medium length or branched |
| **Large** | K, M, Q, H, E | Long or charged sidechains |
| **Bulky** | R, F, Y, W | Aromatic rings or multi-functional groups |

---

## Configuration

### `config.yaml`
```yaml
target_residue: ARG
```

Specify the target amino acid (3-letter code) to analyze. Change to `LYS`, `PHE`, etc. for different residues.

### Environment Variables (Bash)
```bash
export PROTEIN_DIR=protein_structures
export SS_PREDICTOR_BINARY=/usr/local/bin/stride
```

### Command-Line Arguments
```bash
./run_workflow.sh /path/to/proteins /path/to/stride
```

---

## Output Interpretation

### `measurement_angles.tsv`
Tab-separated table with columns:
- `pdb` — PDB identifier
- `left_aa` — single-letter code of preceding residue
- `size_class` — bulk classification (Tiny, Small, Intermediate, Large, Bulky)
- `angle` — signed dihedral angle in degrees (−180 to +180)

### `orientation_density.png`
Density plot showing:
- **X-axis:** Sidechain angle (degrees)
- **Y-axis:** Normalized frequency
- **Colors:** One curve per bulk class
  - Blue (Tiny) → Green (Small) → Orange (Intermediate) → Red (Large) → Dark Red (Bulky)

Interpretation:
- **Shifted distributions:** Bulk class influences average orientation
- **Compressed/expanded distributions:** Bulk class affects orientation flexibility

---

## Troubleshooting

### R/Rscript Not Found
**Error:** `Rscript: command not found` or `not recognized as an internal command`

**Solution:**
- Linux/macOS: Ensure R is installed and in `PATH`: `which Rscript`
- Windows: Script auto-detects R at `C:\Program Files\R\R-*\bin\Rscript.exe`
- Manually specify path in environment or update script

### Missing Python Packages
**Error:** `ModuleNotFoundError: No module named 'pandas'`

**Solution:**
```bash
pip install -r requirements.txt
```

### STRIDE Not Found
**Error:** `FileNotFoundError: [Errno 2] No such file or directory: 'stride'`

**Solution:**
- Install STRIDE from source or package manager
- Linux/macOS: `ln -s /path/to/stride /usr/local/bin/stride`
- Windows: Update `SS_PREDICTOR_BINARY` in PowerShell script

### Snakemake File Not Found
**Error:** `MissingInputException: Missing input files for rule...`

**Solution:** Ensure all `.smk` files are in the project root and run from the project directory.

---

## Advanced Usage

### Analyze Different Residues
Edit `config.yaml`:
```yaml
target_residue: LYS
```

### Parallel Execution
Increase the number of cores:
```bash
./run_workflow.sh --cores 32
```

### Dry Run (No Execution)
Validate workflow without running:
```bash
snakemake -s assign_ss.smk --dry-run
```

### Generate Workflow DAG
Visualize pipeline structure:
```bash
snakemake -s assign_ss.smk --dag | dot -Tpng > dag.png
```

---

## Citation & References

**Pipeline components:**
- Secondary structure: [STRIDE](http://webclu.bio.wzw.tum.de/stride/)
- Workflow: [Snakemake](https://snakemake.readthedocs.io/)
- Structure parsing: [BioPython](https://biopython.org/)
- Visualization: [ggplot2](https://ggplot2.tidyverse.org/)

**Data source:**
- Protein structures: [RCSB Protein Data Bank](https://www.rcsb.org/)

---

## License & Contact

This pipeline is provided as-is for research purposes. Please cite this work if used in publications.

## Output Files

| File | Description |
|------|-------------|
| `results/angles.tsv` | Tab-separated: `pdb`, `left_aa`, `size_class`, `angle` |
| `results/angle_plot.png` | Density curves grouped by size class |

---

## Configuration

To analyse a different residue, update `config.yaml`:

```yaml
target_aa: ARG
```

Any three-letter residue code present in STRIDE ASG output will work.
