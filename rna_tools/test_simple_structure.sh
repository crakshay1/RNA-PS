#!/bin/bash

# Test simple MFE structure drawing

echo "=== Testing Simple MFE Structure Drawing ==="

source .venv/bin/activate

# Get sequence from your best candidate
SEQUENCE=$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)

echo "Using sequence: $SEQUENCE"
echo "Length: ${#SEQUENCE} nucleotides"
echo ""

echo "1. Generating simple MFE structure (arc diagram only)..."
python local_rnafold.py -s "$SEQUENCE" --simple -o simple_structure.png

echo ""
echo "2. Comparing with full dot plot..."
python local_rnafold.py -s "$SEQUENCE" -o full_dotplot.png

echo ""
echo "3. Testing different candidates with simple view..."

for i in {1..3}; do
    CANDIDATE_FILE="designed_sequences/z_tile_tetramer_no_filter_cand${i}.txt"
    if [ -f "$CANDIDATE_FILE" ]; then
        SEQ=$(grep "Sequence:" "$CANDIDATE_FILE" | cut -d' ' -f2)
        MFE=$(grep "MFE:" "$CANDIDATE_FILE" | cut -d' ' -f2)
        echo "Candidate $i (MFE: $MFE)..."
        python local_rnafold.py -s "$SEQ" --simple -o "candidate_${i}_simple.png"
    fi
done

echo ""
echo "=== Generated Files ==="
echo "Simple structure plots:"
echo "- simple_structure.png (clean arc diagram)"
echo "- candidate_1_simple.png, candidate_2_simple.png, candidate_3_simple.png"
echo ""
echo "Full dot plot for comparison:"
echo "- full_dotplot.png (complete analysis with probabilities)"
echo ""
echo "The simple plots show just the MFE structure with:"
echo "- Color-coded arcs for different nesting levels"
echo "- Nucleotide sequence below"
echo "- Position numbers"
echo "- Basic statistics"