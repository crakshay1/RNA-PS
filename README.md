# RNA Structure Processing and Prediction Tool

This directory contains Python scripts for processing RNA structure data, predicting secondary structures, and saving the results.

## `process_rna_data.py`

This script demonstrates how to:
1. Parse both `.cif` and `.pdb` files containing RNA 3D structure data using Biopython.
2. Extract the RNA sequence from the parsed structure.
3. Predict the RNA secondary structure (in dot-bracket notation) and its Minimum Free Energy (MFE) using the ViennaRNA package.
4. Generate arc plot visualizations of the secondary structure using matplotlib.
5. Save the extracted sequence, predicted secondary structure with MFE, and visualization to files.

### Setup

To run this script, you need to set up a Python virtual environment and install the required libraries.

`setup.sh` will create the venv and install required libraries. On linux, `chmod +x setup.sh && ./setup.sh`

### Usage

The script can be run in several ways:

**Default behavior (uses example files):**
```bash
.venv/bin/python3 process_rna_data.py
```

**With specific input file:**
```bash
.venv/bin/python3 process_rna_data.py your_structure.pdb
.venv/bin/python3 process_rna_data.py your_structure.cif
```

**With custom output directory:**
```bash
.venv/bin/python3 process_rna_data.py your_structure.pdb my_results
```

### Output

The script will create an `output` directory (if it doesn't exist) or use the specified output directory. For each processed RNA chain, it will generate three files:

*   `[structure_id]_[chain_id]_sequence.txt`: Contains the extracted RNA sequence.
*   `[structure_id]_[chain_id]_secondary_structure.txt`: Contains the predicted secondary structure in dot-bracket notation and its Minimum Free Energy (MFE).
*   `[structure_id]_[chain_id]_structure_arc_plot.png`: Visualization of the secondary structure as an arc diagram.

## `create_rna_data.py`

This script generates 3D RNA structure models from sequences and secondary structures using external tools. It demonstrates how to:

1. Generate RNA scaffold candidates using the `generate_rna` module
2. Create input files for FARFAR2 (FASTA and secondary structure files)
3. Run FARFAR2 to generate 3D PDB models from sequence and secondary structure constraints
4. Convert PDB files to CIF format using external conversion tools

### Prerequisites

This script requires external tools to be installed:
- **Rosetta/FARFAR2**: For 3D structure prediction
- **pdb2cif** (or equivalent): For PDB to CIF conversion

### Usage

```bash
.venv/bin/python3 create_rna_data.py
```

The script will:
1. Generate 5 candidates for the "z_tile_tetramer" scaffold
2. Save candidates using the `generate_rna` module
3. Use the best candidate to create FARFAR2 input files
4. Run FARFAR2 to generate 3 PDB structure models
5. Convert the first PDB to CIF format

### Output

The script creates several directories and files:
- `farfar_inputs/`: Contains FASTA and secondary structure files for FARFAR2
- `farfar_outputs/`: Contains generated PDB files from FARFAR2
- Generated CIF files converted from PDB models

### Functions

- `write_farfar2_inputs()`: Creates FASTA and secondary structure files for FARFAR2
- `run_farfar2()`: Executes FARFAR2 to generate 3D models
- `convert_pdb_to_cif()`: Converts PDB files to CIF format
- `write_rnacomposer_input()`: Writes input for RNAComposer (alternative tool)

## Pipeline Flow

generate_rna.py → 2. create_rna_data.py → 3. process_rna_data.py

- Step 1: RNA Design (generate_rna.py)
  - Generates RNA scaffold candidates with sequences and secondary structures
  - Creates multiple candidate designs for a given scaffold
  - Outputs: RNA sequences + predicted/target secondary structures

- Step 2: 3D Structure Generation (create_rna_data.py)
  - Takes the best candidate from generate_rna.py
  - Calls gr.generate_candidates_for_scaffold() and gr.save_candidates()
  - Converts sequence + secondary structure → 3D PDB models via FARFAR2
  - Converts PDB → CIF format
  - Outputs: 3D structure files (.pdb, .cif)

- Step 3: Analysis & Visualization (process_rna_data.py)
  - Takes the generated .pdb/.cif files from Step 2
  - Parses 3D structures to extract sequences
  - Predicts secondary structures using ViennaRNA
  - Generates arc plot visualizations
  - Outputs: Sequence files, secondary structure predictions, visualizations

Data Flow:

RNA Scaffold Concept
    ↓ (generate_rna.py)
Sequence + Secondary Structure Candidates
    ↓ (create_rna_data.py)
3D Structure Models (.pdb/.cif)
    ↓ (process_rna_data.py)
Analysis Results + Visualizations



## References

This script utilizes several open-source libraries and a public dataset. Please refer to their documentation and publications for more details.

### Biopython

*   **Documentation:** [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
*   **Publication:** Cock, P. J. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.

### ViennaRNA

*   **Documentation:** [ViennaRNA Python API](https://pypi.org/project/ViennaRNA/)
*   **Publication:** Lorenz, R., et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6(1), 26.

