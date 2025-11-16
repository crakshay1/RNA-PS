#!/bin/bash


echo "=== Testing Dot Plot Visualization (like RNAfold web interface) ==="

source .venv/bin/activate

SEQUENCE=$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)

echo "Using sequence: $SEQUENCE"
echo "Length: ${#SEQUENCE} nucleotides"
echo ""

echo "1. Generating dot plot with base pair probabilities..."
python local_rnafold.py -s "$SEQUENCE" -o dotplot_result.png

echo ""
echo "2. Testing with different temperature (25°C)..."
python local_rnafold.py -s "$SEQUENCE" -t 25 -o dotplot_25C.png

echo ""
echo "3. Testing with higher temperature (60°C)..."
python local_rnafold.py -s "$SEQUENCE" -t 60 -o dotplot_60C.png
