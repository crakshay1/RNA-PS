# RNA Structure Processing and Prediction Tool

This directory contains Python scripts for processing RNA structure data from CIF files, predicting secondary structures, and saving the results.

## `process_rna_data.py`

This script demonstrates how to:
1. Parse a `.cif` file containing RNA 3D structure data using Biopython.
2. Extract the RNA sequence from the parsed structure.
3. Predict the RNA secondary structure (in dot-bracket notation) and its Minimum Free Energy (MFE) using the ViennaRNA package.
4. Save the extracted sequence and the predicted secondary structure with MFE to text files.

### Setup

To run this script, you need to set up a Python virtual environment and install the required libraries.

1.  **Create a virtual environment (if you haven't already):**
    ```bash
    python3 -m venv .venv
    ```

2.  **Activate the virtual environment:**
    ```bash
    source .venv/bin/activate
    ```
    (On Windows, use `.venv\Scripts\activate`)

3.  **Install the required Python packages:**
    ```bash
    .venv/bin/pip install biopython viennarna
    ```

### Usage

The script is designed to process a single `.cif` file. You can modify the `example_cif_file` variable in the `if __name__ == "__main__":` block to point to any `.cif` file in your `rna3db-mmcifs` database.
It is designed to work with [rna3db](https://github.com/marcellszi/rna3db/releases/tag/2025-10-01-incremental-release)

**Example:**

```bash
.venv/bin/python3 code/rna_tools/process_rna_data.py
```

### Database Structure

The script expects the `rna3db-mmcifs` database to be located at `code/databases/rna3db-mmcifs` relative to the project root. The database has a nested structure:

```
code/databases/rna3db-mmcifs/
├── train_set/
│   ├── component_1/
│   │   ├── 1c2w_B/
│   │   │   └── 1c2w_B.cif
│   │   ├── ...
│   ├── component_X/
│   │   └── ...
└── test_set/
    └── ...
```

To process a different file, update the `example_cif_file` variable in `process_rna_data.py` accordingly.

### Output

The script will create an `output` directory (if it doesn't exist) in the project root. For each processed RNA chain, it will generate two text files:

*   `[structure_id]_[chain_id]_sequence.txt`: Contains the extracted RNA sequence.
*   `[structure_id]_[chain_id]_secondary_structure.txt`: Contains the predicted secondary structure in dot-bracket notation and its Minimum Free Energy (MFE).


## References

This script utilizes several open-source libraries and a public dataset. Please refer to their documentation and publications for more details.

### Biopython

*   **Documentation:** [Biopython Tutorial and Cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
*   **Publication:** Cock, P. J. A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.

### ViennaRNA

*   **Documentation:** [ViennaRNA Python API](https://pypi.org/project/ViennaRNA/)
*   **Publication:** Lorenz, R., et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6(1), 26.

### rna3db-mmcifs Dataset

*   **Publication:** Szikszai, M., et al. (2024). RNA3DB: A structurally-dissimilar dataset split for training and benchmarking deep learning models for RNA structure prediction. *Journal of Molecular Biology*, 436(17).
