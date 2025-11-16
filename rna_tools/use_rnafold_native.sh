#!/bin/bash

echo "=== Using Native RNAfold Command-Line Tool ==="

# Get sequence from your best candidate
SEQUENCE=$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)
MFE_EXPECTED=$(grep "MFE:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)

echo "Using sequence: $SEQUENCE"
echo "Expected MFE: $MFE_EXPECTED"
echo ""

# Create a temporary FASTA file
echo ">z_tile_tetramer_candidate1" > temp_sequence.fa
echo "$SEQUENCE" >> temp_sequence.fa

echo "1. Basic MFE calculation..."
echo "Command: RNAfold < temp_sequence.fa"
RNAfold < temp_sequence.fa
echo ""

echo "2. MFE + Partition function + Base pair probabilities..."
echo "Command: RNAfold -p < temp_sequence.fa"
RNAfold -p < temp_sequence.fa
echo ""

echo "3. Generate plots with specific temperature (25°C)..."
echo "Command: RNAfold -p -T 25 < temp_sequence.fa"
RNAfold -p -T 25 < temp_sequence.fa
echo ""

echo "4. Generate plots with constraints (if you want to force some base pairs)..."
echo "Example with no specific constraints:"
RNAfold -p -T 37 < temp_sequence.fa
echo ""

echo "5. Testing multiple sequences..."
echo "Creating multi-sequence file..."

# Create multi-sequence FASTA
cat > multi_sequences.fa << EOF
>candidate1
$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)
>candidate2
$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand2.txt | cut -d' ' -f2)
>candidate3
$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand3.txt | cut -d' ' -f2)
EOF

echo "Processing multiple sequences..."
RNAfold -p < multi_sequences.fa

echo ""
echo "=== Generated Files ==="
echo "PostScript files created:"
ls -la *.ps 2>/dev/null || echo "No .ps files found"
echo ""

echo "=== Converting PostScript to PNG (if you prefer) ==="
echo "To convert PS to PNG, install ghostscript and use:"
echo "  sudo apt install ghostscript"
echo "  gs -sDEVICE=png16m -sOutputFile=structure.png -r300 -dNOPAUSE -dBATCH rna_ss.ps"
echo "  gs -sDEVICE=png16m -sOutputFile=dotplot.png -r300 -dNOPAUSE -dBATCH rna_dp.ps"
echo ""

echo "=== Key RNAfold Options ==="
echo "-p          : Calculate partition function and base pair probabilities"
echo "-T <temp>   : Set temperature (default 37°C)"
echo "-d2         : Allow dangles on both sides of helices"
echo "--noPS      : Skip PostScript structure plot"
echo "--noDP      : Skip PostScript dot plot"
echo "-C          : Read constraints from input"
echo "--MEA       : Calculate maximum expected accuracy structure"

# Cleanup
rm -f temp_sequence.fa