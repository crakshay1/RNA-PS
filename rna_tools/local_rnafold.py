#!/usr/bin/env python3

import RNA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
import sys
import os
from datetime import datetime

def fold_rna_sequence(sequence, temperature=37.0, constraints=None):
    """
    Fold RNA sequence using ViennaRNA with optional constraints.

    Args:
        sequence: RNA sequence string
        temperature: Temperature in Celsius (default 37.0)
        constraints: Optional structure constraints in dot-bracket notation

    Returns:
        tuple: (secondary_structure, mfe_energy, ensemble_energy, centroid_structure)
    """
    # Set temperature
    RNA.cvar.temperature = temperature

    # Create fold compound
    fc = RNA.fold_compound(sequence)

    # Apply constraints if provided
    if constraints:
        fc.constraints_add(constraints, RNA.CONSTRAINT_DB_DEFAULT)

    # Calculate MFE structure
    mfe_structure, mfe_energy = fc.mfe()

    # Calculate partition function and ensemble properties
    ensemble_energy = fc.pf()[1]

    # Get centroid structure
    centroid_structure, centroid_distance = fc.centroid()

    # Calculate base pair probabilities
    bp_prob_matrix = fc.bpp()

    return {
        'sequence': sequence,
        'mfe_structure': mfe_structure,
        'mfe_energy': mfe_energy,
        'ensemble_energy': ensemble_energy,
        'centroid_structure': centroid_structure,
        'centroid_distance': centroid_distance,
        'bp_probabilities': bp_prob_matrix,
        'temperature': temperature
    }

def plot_dot_plot(result, output_path=None, show_plot=True):
    """
    Create dot plot visualization showing base pair probabilities like RNAfold web interface.
    """
    sequence = result['sequence']
    bp_matrix = result['bp_probabilities']
    structure = result['mfe_structure']
    n = len(sequence)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Base pair probability dot plot (upper triangle)
    ax1.set_title(f"Base Pair Probabilities - Length: {n}")

    # Create probability matrix for visualization
    prob_matrix = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            prob = bp_matrix[i][j]
            if prob > 0.01:  # Only show significant probabilities
                prob_matrix[i-1][j-1] = prob

    # Custom colormap: white -> light blue -> dark blue -> red (like web interface)
    colors = ['white', '#e6f3ff', '#99d6ff', '#4da6ff', '#0066cc', '#ff3300']
    cmap = mcolors.LinearSegmentedColormap.from_list('probability', colors, N=100)

    # Plot upper triangle only
    mask = np.triu(np.ones_like(prob_matrix), k=1)
    prob_matrix_masked = np.where(mask, prob_matrix, np.nan)

    im1 = ax1.imshow(prob_matrix_masked, cmap=cmap, vmin=0, vmax=1, origin='upper')
    ax1.set_xlabel('Position j')
    ax1.set_ylabel('Position i')

    # Add grid and ticks
    ax1.set_xticks(range(0, n, 10))
    ax1.set_yticks(range(0, n, 10))
    ax1.grid(True, alpha=0.3)

    # Add colorbar
    cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
    cbar1.set_label('Base Pair Probability')

    # Plot 2: MFE structure as dot plot (showing actual base pairs)
    ax2.set_title(f"MFE Structure - Energy: {result['mfe_energy']:.2f} kcal/mol")

    # Create matrix for MFE structure
    mfe_matrix = np.zeros((n, n))
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            mfe_matrix[i][j] = 1.0
            mfe_matrix[j][i] = 1.0

    # Plot MFE structure
    mask = np.triu(np.ones_like(mfe_matrix), k=1)
    mfe_matrix_masked = np.where(mask, mfe_matrix, np.nan)

    im2 = ax2.imshow(mfe_matrix_masked, cmap='Reds', vmin=0, vmax=1, origin='upper')
    ax2.set_xlabel('Position j')
    ax2.set_ylabel('Position i')
    ax2.set_xticks(range(0, n, 10))
    ax2.set_yticks(range(0, n, 10))
    ax2.grid(True, alpha=0.3)

    # Plot 3: Arc diagram
    ax3.set_title("Arc Diagram Representation")
    positions = range(n)
    ax3.plot(positions, [0] * n, 'ko-', markersize=3, linewidth=1)

    # Add nucleotide labels
    for i, nucleotide in enumerate(sequence):
        ax3.text(i, -0.15, nucleotide, ha='center', va='top', fontsize=6)
        if i % 10 == 0:  # Position numbers every 10
            ax3.text(i, 0.15, str(i+1), ha='center', va='bottom', fontsize=6, color='gray')

    # Draw base pairs as arcs
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            # Draw arc
            center = (i + j) / 2
            width = i - j
            height = width * 0.3
            arc = patches.Arc((center, 0), width=width, height=height,
                            theta1=0, theta2=180, edgecolor='blue', linewidth=1.5)
            ax3.add_patch(arc)

    ax3.set_xlim(-1, n)
    ax3.set_ylim(-0.5, n/4)
    ax3.set_xlabel('Position')
    ax3.grid(True, alpha=0.3)

    # Plot 4: Statistics and information
    ax4.axis('off')
    info_text = f"""
Sequence: {sequence[:50]}{"..." if len(sequence) > 50 else ""}
Length: {n} nucleotides

MFE Structure: {structure[:50]}{"..." if len(structure) > 50 else ""}
MFE Energy: {result['mfe_energy']:.2f} kcal/mol

Ensemble Energy: {result['ensemble_energy']:.2f} kcal/mol
Centroid Structure: {result['centroid_structure'][:50]}{"..." if len(result['centroid_structure']) > 50 else ""}
Centroid Distance: {result['centroid_distance']:.2f}

Temperature: {result['temperature']}°C
Algorithm: ViennaRNA Package
Computed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Base Pairs in MFE: {structure.count('(')}
GC Content: {((sequence.count('G') + sequence.count('C')) / len(sequence) * 100):.1f}%
    """

    ax4.text(0.05, 0.95, info_text, transform=ax4.transAxes, fontsize=9,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Dot plot saved to {output_path}")

    if show_plot:
        plt.show()

    return fig

    # Plot 1: Arc diagram (similar to web interface)
    ax1.set_title(f"RNA Secondary Structure - MFE: {result['mfe_energy']:.2f} kcal/mol")

    # Draw sequence positions
    positions = range(len(sequence))
    ax1.plot(positions, [0] * len(sequence), 'ko-', markersize=4, linewidth=1)

    # Add nucleotide labels
    for i, nucleotide in enumerate(sequence):
        ax1.text(i, -0.1, nucleotide, ha='center', va='top', fontsize=8)
        ax1.text(i, 0.1, str(i+1), ha='center', va='bottom', fontsize=6, color='gray')

    # Draw base pairs as arcs
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            # Draw arc
            center = (i + j) / 2
            radius = (i - j) / 2
            height = radius * 0.3
            arc = patches.Arc((center, 0), width=i-j, height=height,
                            theta1=0, theta2=180, edgecolor='blue', linewidth=1.5)
            ax1.add_patch(arc)

    ax1.set_xlim(-1, len(sequence))
    ax1.set_ylim(-0.5, len(sequence)/4)
    ax1.set_xlabel('Position')
    ax1.set_ylabel('')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Energy and statistics
    ax2.axis('off')
    info_text = f"""
Sequence Length: {len(sequence)} nucleotides
Sequence: {sequence}

MFE Structure: {structure}
MFE Energy: {result['mfe_energy']:.2f} kcal/mol

Ensemble Energy: {result['ensemble_energy']:.2f} kcal/mol
Centroid Structure: {result['centroid_structure']}
Centroid Distance: {result['centroid_distance']:.2f}

Temperature: {result['temperature']}°C
Folding Algorithm: ViennaRNA {RNA.__version__ if hasattr(RNA, '__version__') else 'N/A'}
Computed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    """

    ax2.text(0.05, 0.95, info_text, transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_path}")

    if show_plot:
        plt.show()

    return fig

def plot_simple_structure(result, output_path=None, show_plot=True):
    """
    Create a simple, clean MFE structure visualization - just the arc diagram.
    """
    sequence = result['sequence']
    structure = result['mfe_structure']
    n = len(sequence)

    fig, ax = plt.subplots(1, 1, figsize=(max(12, n/5), 6))

    ax.set_title(f"RNA Secondary Structure\nMFE: {result['mfe_energy']:.2f} kcal/mol",
                fontsize=14, pad=20)

    # Draw sequence positions as line
    positions = range(n)
    ax.plot(positions, [0] * n, 'ko-', markersize=4, linewidth=1, color='black')

    # Add nucleotide labels below the line
    for i, nucleotide in enumerate(sequence):
        ax.text(i, -0.2, nucleotide, ha='center', va='top', fontsize=10,
               fontweight='bold', color='darkblue')

    # Add position numbers above every 10th position
    for i in range(0, n, 10):
        ax.text(i, 0.3, str(i+1), ha='center', va='bottom', fontsize=8,
               color='gray', alpha=0.8)

    # Draw base pairs as arcs
    stack = []
    colors = ['blue', 'red', 'green', 'orange', 'purple']  # Different colors for nested levels
    level = 0

    for i, char in enumerate(structure):
        if char == '(':
            stack.append((i, level))
            level += 1
        elif char == ')' and stack:
            j, pair_level = stack.pop()
            level -= 1

            # Draw arc
            center = (i + j) / 2
            width = i - j
            height = width * 0.4 + pair_level * 0.5  # Adjust height based on nesting level

            arc_color = colors[pair_level % len(colors)]
            arc = patches.Arc((center, 0), width=width, height=height,
                            theta1=0, theta2=180, edgecolor=arc_color,
                            linewidth=2, alpha=0.8)
            ax.add_patch(arc)

    # Set limits and styling
    ax.set_xlim(-2, n + 1)
    max_height = max(10, n/4)
    ax.set_ylim(-1, max_height)

    # Remove y-axis and ticks
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Style x-axis
    ax.set_xlabel('Nucleotide Position', fontsize=12)
    ax.set_xticks(range(0, n, max(1, n//20)))  # Show reasonable number of ticks

    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='x')

    # Add statistics as text
    stats_text = f"Length: {n} nt  |  Base pairs: {structure.count('(')}  |  " \
                f"GC content: {((sequence.count('G') + sequence.count('C')) / n * 100):.1f}%  |  " \
                f"Temperature: {result['temperature']}°C"

    ax.text(0.5, 0.02, stats_text, transform=ax.transAxes, ha='center',
           fontsize=10, bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Simple structure plot saved to {output_path}")

    if show_plot:
        plt.show()

    return fig

def calculate_base_pair_probabilities(result, min_prob=0.01):
    """
    Extract and display base pair probabilities above threshold.
    """
    bp_probs = []
    bp_matrix = result['bp_probabilities']
    sequence_length = len(result['sequence'])

    for i in range(1, sequence_length + 1):
        for j in range(i + 1, sequence_length + 1):
            prob = bp_matrix[i][j]
            if prob > min_prob:
                bp_probs.append((i-1, j-1, prob))  # Convert to 0-indexed

    return sorted(bp_probs, key=lambda x: x[2], reverse=True)

def main():
    parser = argparse.ArgumentParser(
        description="Local RNAfold - RNA secondary structure prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python local_rnafold.py -s GGGAAACCC
  python local_rnafold.py -s GGGAAACCC -t 25 -o result.png
  python local_rnafold.py -f sequence.txt --constraints "(((...)))"
        """
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--sequence', help='RNA sequence to fold')
    group.add_argument('-f', '--file', help='File containing RNA sequence')

    parser.add_argument('-t', '--temperature', type=float, default=37.0,
                       help='Temperature in Celsius (default: 37.0)')
    parser.add_argument('-c', '--constraints',
                       help='Structure constraints in dot-bracket notation')
    parser.add_argument('-o', '--output',
                       help='Output file for plot (PNG format)')
    parser.add_argument('--no-plot', action='store_true',
                       help='Skip plotting, only show results')
    parser.add_argument('--simple', action='store_true',
                       help='Generate simple arc diagram instead of full dot plot')
    parser.add_argument('--min-prob', type=float, default=0.01,
                       help='Minimum base pair probability to display (default: 0.01)')

    args = parser.parse_args()

    # Get sequence
    if args.sequence:
        sequence = args.sequence.upper()
    else:
        with open(args.file, 'r') as f:
            sequence = f.read().strip().replace('\n', '').replace(' ', '').upper()

    # Validate sequence
    valid_nucleotides = set('AUGC')
    if not all(nuc in valid_nucleotides for nuc in sequence):
        print("Error: Invalid nucleotides in sequence. Only A, U, G, C allowed.")
        sys.exit(1)

    print(f"Folding sequence: {sequence}")
    print(f"Length: {len(sequence)} nucleotides")
    print(f"Temperature: {args.temperature}°C")

    # Fold the sequence
    result = fold_rna_sequence(sequence, args.temperature, args.constraints)

    # Display results
    print("\n=== RESULTS ===")
    print(f"MFE Structure: {result['mfe_structure']}")
    print(f"MFE Energy: {result['mfe_energy']:.2f} kcal/mol")
    print(f"Ensemble Energy: {result['ensemble_energy']:.2f} kcal/mol")
    print(f"Centroid Structure: {result['centroid_structure']}")
    print(f"Centroid Distance: {result['centroid_distance']:.2f}")

    # Show base pair probabilities
    bp_probs = calculate_base_pair_probabilities(result, args.min_prob)
    if bp_probs:
        print(f"\n=== BASE PAIR PROBABILITIES (> {args.min_prob}) ===")
        for i, j, prob in bp_probs[:20]:  # Show top 20
            print(f"({i+1:3d}, {j+1:3d}): {prob:.4f}")
        if len(bp_probs) > 20:
            print(f"... and {len(bp_probs) - 20} more pairs")

    # Create plot
    if not args.no_plot:
        if args.simple:
            plot_simple_structure(result, args.output, show_plot=not args.output)
        else:
            plot_dot_plot(result, args.output, show_plot=not args.output)

if __name__ == "__main__":
    main()