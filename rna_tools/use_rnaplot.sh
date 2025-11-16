#!/bin/bash

echo "=== Using RNAplot for RNA Structure Visualization ==="

# Get sequence and structure from RNAfold output
SEQUENCE=$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand1.txt | cut -d' ' -f2)

echo "Using sequence: $SEQUENCE"
echo ""

echo "1. Generate structure with RNAfold and plot with RNAplot..."
echo "$SEQUENCE" | RNAfold | RNAplot

echo "2. Create custom named plots..."
echo -e ">my_rna_structure" > temp_input.txt
echo "$SEQUENCE" | RNAfold >> temp_input.txt
RNAplot < temp_input.txt

echo "3. Generate different plot formats..."
# RNAplot with different output formats
echo "$SEQUENCE" | RNAfold | RNAplot --format=svg -o structure.svg
echo "$SEQUENCE" | RNAfold | RNAplot --format=png -o structure.png
echo "$SEQUENCE" | RNAfold | RNAplot --format=pdf -o structure.pdf

echo "4. Advanced plotting options..."
# With sequence annotations
echo "$SEQUENCE" | RNAfold | RNAplot --pre="set title 'Z-Tile Tetramer RNA Structure'"

echo "5. Batch process multiple candidates..."
for i in {1..3}; do
    CANDIDATE_FILE="designed_sequences/z_tile_tetramer_no_filter_cand${i}.txt"
    if [ -f "$CANDIDATE_FILE" ]; then
        SEQ=$(grep "Sequence:" "$CANDIDATE_FILE" | cut -d' ' -f2)
        MFE=$(grep "MFE:" "$CANDIDATE_FILE" | cut -d' ' -f2)

        echo -e ">candidate_${i}_MFE_${MFE}" > candidate_${i}_input.txt
        echo "$SEQ" | RNAfold >> candidate_${i}_input.txt
        RNAplot < candidate_${i}_input.txt

        # Convert to PNG
        gs -sDEVICE=png16m -sOutputFile=candidate_${i}_structure.png -r300 -dNOPAUSE -dBATCH candidate_${i}_MFE_${MFE}_ss.ps
    fi
done

echo ""
echo "=== Generated Files ==="
ls -la *.ps *.png *.svg *.pdf 2>/dev/null | grep -E "(structure|candidate.*ss)"

echo ""
echo "=== RNAplot Commands Summary ==="
echo "Basic usage:"
echo "  echo \"SEQUENCE\" | RNAfold | RNAplot"
echo ""
echo "Different formats:"
echo "  RNAplot --format=svg -o output.svg"
echo "  RNAplot --format=png -o output.png"
echo "  RNAplot --format=pdf -o output.pdf"
echo ""
echo "View PostScript files:"
echo "  evince file.ps"
echo "  gv file.ps"
echo "  okular file.ps"

# Cleanup
rm -f temp_input.txt candidate_*_input.txt