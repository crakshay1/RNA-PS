# RNA Structure Processing - PDB and CIF Support

The `process_rna_data.py` script now supports both PDB and CIF file formats for RNA structure processing.

## Usage Examples

### Command Line Usage

```bash
# Process a PDB file
python process_rna_data.py my_rna_structure.pdb

# Process a CIF file (original functionality)
python process_rna_data.py my_rna_structure.cif

# Specify custom output directory
python process_rna_data.py my_rna_structure.pdb custom_output_dir

# Using with virtual environment
source .venv/bin/activate && python process_rna_data.py structure.pdb
```

### Python API Usage

```python
import process_rna_data as prd

# Unified function - automatically detects file type
prd.process_structure_file("structure.pdb", "output_dir")
prd.process_structure_file("structure.cif", "output_dir")

# Specific functions
prd.process_pdb_file("structure.pdb", "output_dir")
prd.process_cif_file("structure.cif", "output_dir")  # Legacy
```

## What the Script Does

1. **Parses structure files** using appropriate BioPython parser (PDB or CIF)
2. **Extracts RNA sequences** from all chains
3. **Predicts secondary structure** using ViennaRNA
4. **Calculates MFE** (Minimum Free Energy)
5. **Saves results** as text files with sequence and structure info
6. **Generates visualizations** as PNG arc diagrams

## Output Files

For each RNA chain found, the script creates:

- `{structure_id}_{chain_id}_sequence.txt` - RNA sequence
- `{structure_id}_{chain_id}_secondary_structure.txt` - Structure + MFE
- `{structure_id}_{chain_id}_structure_arc_plot.png` - Visualization

## Supported File Formats

| Format | Extension | Parser Used | Notes |
|--------|-----------|-------------|-------|
| PDB | `.pdb` | PDBParser | Standard Protein Data Bank format |
| CIF | `.cif` | MMCIFParser | Crystallographic Information Format |

## Integration with Your Workflow

### Using with Generated Sequences

```bash
# Convert your generated sequences to 3D models (PDB format)
# Then process them back to analyze the structures

# 1. Generate sequences with MFE filtering
python generate_rna.py

# 2. Create 3D models (using FARFAR2 or similar)
python create_rna_data.py

# 3. Process the resulting PDB files
python process_rna_data.py output_model.pdb
```

### Batch Processing Multiple Files

```bash
# Process all PDB files in a directory
for pdb_file in *.pdb; do
    echo "Processing $pdb_file..."
    python process_rna_data.py "$pdb_file" "results_$(basename "$pdb_file" .pdb)"
done

# Process all CIF files
for cif_file in *.cif; do
    echo "Processing $cif_file..."
    python process_rna_data.py "$cif_file" "results_$(basename "$cif_file" .cif)"
done
```

## Error Handling

The script handles common issues:

- **File not found**: Clear error message
- **Unsupported format**: Only .pdb and .cif are accepted
- **Parsing errors**: BioPython exceptions are caught
- **No RNA found**: Skips chains without standard RNA nucleotides (A, U, G, C)
- **Visualization errors**: Continues processing even if plotting fails

## Testing

```bash
# Test with sample files
python test_pdb_processing.py

# Test the modified script
python process_rna_data.py test_structure.pdb test_output
```

This enhancement makes your pipeline compatible with both the original CIF database and any PDB files you generate from your RNA design workflow!