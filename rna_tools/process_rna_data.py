
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.PDB import MMCIFParser
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

def process_cif_file(cif_file_path, output_dir="output"):
    """
    Parses a .cif file, predicts secondary structure, saves the data, and generates a visualization.
    """
    if not os.path.exists(cif_file_path):
        print(f"Error: File not found at {cif_file_path}")
        return

    parser = MMCIFParser()
    try:
        structure = parser.get_structure("RNA_structure", cif_file_path)
    except PDBException as e:
        print(f"Error parsing {cif_file_path}: {e}")
        return

    print(f"Successfully parsed {cif_file_path}")
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


if __name__ == "__main__":

    example_cif_file = "code/databases/rna3db-mmcifs/train_set/component_1/1c2w_B/1c2w_B.cif"
    
    try:
        import Bio
        import RNA
        import matplotlib
    except ImportError:
        print("Required libraries (biopython, viennarna, matplotlib) are not installed. Please install them using '.venv/bin/pip install biopython viennarna matplotlib'")
    else:
        process_cif_file(example_cif_file)
