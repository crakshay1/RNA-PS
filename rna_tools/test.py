#!/usr/bin/env python

# Make sure you have installed ViennaRNA via conda
import RNA
import os

# --- Your RNA sequence ---
sequence = "GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUUGUACGGGAUGUCCGGCGUUCGAGUCCGGCAGCUCCCA"

print(f"Sequence: {sequence}\n")

# --- 1. Perform All Calculations ---

# Create a "fold compound" object for the sequence
fc = RNA.fold_compound(sequence)

# Calculate the MFE (Minimum Free Energy) structure
(mfe_struct, mfe) = fc.mfe()
print(f"MFE Structure:    {mfe_struct}")

# Calculate the partition function (PF)
# This is ESSENTIAL for probabilities and the centroid
(pf_struct, pf_energy) = fc.pf()

# Get the Centroid structure
(centroid_struct, distance) = fc.centroid()
print(f"Centroid Structure: {centroid_struct}\n")

# --- Get the Base-Pair Probability Matrix ---
# This is the data used for the color-coding.
# We call fc.bpp() *after* fc.pf()
bppm = fc.bpp() 

# --- 2. Plot the MFE Structure ---
# We pass the probability matrix 'bppm' using the 'bpp=' keyword.
mfe_filename = "mfe_plot.svg"
RNA.svg_rna_plot(sequence, mfe_struct, mfe_filename, bpp=bppm) # <-- CORRECTED
print(f"Successfully created MFE plot: {os.path.abspath(mfe_filename)}")


# --- 3. Plot the Centroid Structure ---
# We pass the *same* probability matrix 'bppm'
centroid_filename = "centroid_plot.svg"
RNA.svg_rna_plot(sequence, centroid_struct, centroid_filename, bpp=bppm) # <-- CORRECTED
print(f"Successfully created Centroid plot: {os.path.abspath(centroid_filename)}")