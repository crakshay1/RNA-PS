
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBException
import RNA 

def plot_arc_diagram(ss, sequence, output_png_path):
    """
    Generates and saves an arc plot of the RNA secondary structure using matplotlib.
    """
    fig, ax = plt.subplots()
    ax.plot(range(len(sequence)), [0] * len(sequence), 'o-', color='gray', markersize=2)

    stack = []
    for i, char in enumerate(ss):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                arc = patches.Arc(((i + j) / 2, 0), width=i - j, height=(i - j) / 2,
                              theta1=0, theta2=180, edgecolor='blue', lw=1)
                ax.add_patch(arc)

    ax.set_yticks([])
    ax.set_xticks(range(len(sequence)))
    ax.set_xticklabels(list(sequence))
    ax.set_xlim(-1, len(sequence))
    ax.set_ylim(-5, len(sequence) / 2)
    ax.set_aspect('equal')
    plt.title("RNA Secondary Structure Arc Plot")
    plt.savefig(output_png_path)
    plt.close()

def process_structure_file(file_path, output_dir="output"):
    """
    Parses a .cif or .pdb file, predicts secondary structure, saves the data, and generates a visualization.

    Args:
        file_path (str): Path to PDB or CIF structure file
        output_dir (str): Directory to save output files (default: "output")

    Example usage:
        process_structure_file("Predict.pdb")
        process_structure_file("structure.cif", "my_results")
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    # Determine file type and choose appropriate parser
    file_extension = file_path.lower().split('.')[-1]

    if file_extension == 'cif':
        parser = MMCIFParser()
        print(f"Using MMCIFParser for {file_path}")
    elif file_extension == 'pdb':
        parser = PDBParser()
        print(f"Using PDBParser for {file_path}")
    else:
        print(f"Error: Unsupported file format '{file_extension}'. Only .cif and .pdb files are supported.")
        return

    try:
        structure = parser.get_structure("RNA_structure", file_path)
    except PDBException as e:
        print(f"Error parsing {file_path}: {e}")
        return

    print(f"Successfully parsed {file_path}")
    print(f"Structure ID: {structure.id}")

    os.makedirs(output_dir, exist_ok=True)

    for chain in structure.get_chains():
        print(f"Processing chain {chain.id}...")

        # Extract sequence
        sequence = ""
        for residue in chain.get_residues():
            if residue.get_resname().strip() in ("A", "U", "G", "C"):
                sequence += residue.get_resname().strip()

        if not sequence:
            print(f"No standard RNA sequence found for chain {chain.id}")
            continue
        
        """
        if len(sequence) > 300:
            print("Sequence is too long for arc plot visualization, skipping.")
            continue
        """
        print(f"Sequence: {sequence}")

        (ss, energy) = RNA.fold(sequence)
        print(f"Secondary Structure: {ss} (MFE: {energy:.2f})")

        sequence_output_path = os.path.join(output_dir, f"{structure.id}_{chain.id}_sequence.txt")
        with open(sequence_output_path, "w") as f:
            f.write(sequence)
        print(f"Sequence saved to {sequence_output_path}")

        structure_output_path = os.path.join(output_dir, f"{structure.id}_{chain.id}_secondary_structure.txt")
        with open(structure_output_path, "w") as f:
            f.write(f"Sequence: {sequence}\n")
            f.write(f"Secondary Structure: {ss}\n")
            f.write(f"MFE: {energy:.2f}\n")
        print(f"Secondary structure and MFE saved to {structure_output_path}")

        print("Generating visualization...")
        try:
            png_output_path = os.path.join(output_dir, f"{structure.id}_{chain.id}_structure_arc_plot.png")
            plot_arc_diagram(ss, sequence, png_output_path)
            print(f"Visualization saved to {png_output_path}")
        except Exception as e:
            print(f"An unexpected error occurred during visualization: {e}")

# Legacy function for backward compatibility
def process_cif_file(cif_file_path, output_dir="output"):
    """
    Legacy function - use process_structure_file() instead.
    Parses a .cif file, predicts secondary structure, saves the data, and generates a visualization.
    """
    return process_structure_file(cif_file_path, output_dir)

def process_pdb_file(pdb_file_path, output_dir="output2"):
    """
    Parses a .pdb file, predicts secondary structure, saves the data, and generates a visualization.
    """
    return process_structure_file(pdb_file_path, output_dir)

if __name__ == "__main__":

    import sys

    # Check for command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "output"
    else:
        # Default example files - try PDB first, then CIF as fallback
        example_pdb_file = "Predict.pdb"
        example_cif_file = "/home/paul/Documents/EurotechBioHack/code/databases/rna3db-mmcifs/train_set/component_1/1c2w_B/1c2w_B.cif"

        # Use PDB file if it exists, otherwise fall back to CIF
        if os.path.exists(example_pdb_file):
            input_file = example_pdb_file
            print(f"Using default PDB file: {example_pdb_file}")
        elif os.path.exists(example_cif_file):
            input_file = example_cif_file
            print(f"PDB file not found, using default CIF file: {example_cif_file}")
        else:
            print("Error: No default example files found.")
            print(f"Please provide a PDB or CIF file as argument:")
            print(f"  python {sys.argv[0]} your_structure.pdb")
            print(f"  python {sys.argv[0]} your_structure.cif")
            sys.exit(1)

        output_dir = "output"

    try:
        import Bio
        import RNA
        import matplotlib
    except ImportError:
        print("Required libraries (biopython, viennarna, matplotlib) are not installed.")
        print("Please install them using: .venv/bin/pip install biopython viennarna matplotlib")
    else:
        print(f"Processing file: {input_file}")
        print(f"Output directory: {output_dir}")
        print("=" * 50)

        # Use the unified function that handles both PDB and CIF
        process_structure_file(input_file, output_dir)
