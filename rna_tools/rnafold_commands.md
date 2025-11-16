# RNAfold Native Command-Line Usage

## Quick Commands with Your Sequences

### 1. Basic MFE Structure
```bash
# Using your best sequence directly
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" | RNAfold
```

### 2. MFE + Dot Plot (like web interface)
```bash
# Generate structure + base pair probability plots
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" | RNAfold -p
```

### 3. With FASTA Format
```bash
# Create FASTA file
echo ">my_rna_sequence" > my_sequence.fa
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" >> my_sequence.fa

# Fold it
RNAfold -p < my_sequence.fa
```

### 4. Different Temperatures
```bash
# Room temperature (25°C)
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" | RNAfold -p -T 25

# Body temperature (37°C - default)
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" | RNAfold -p -T 37

# High temperature (60°C)
echo "CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC" | RNAfold -p -T 60
```

## Output Files Generated

**With `-p` flag, RNAfold creates:**
- `rna_ss.ps` - Secondary structure plot (PostScript)
- `rna_dp.ps` - Dot plot with base pair probabilities (PostScript)

**For named sequences (FASTA), files are named:**
- `<sequence_name>_ss.ps` - Structure plot
- `<sequence_name>_dp.ps` - Dot plot

## Convert PostScript to PNG/PDF

```bash
# Install ghostscript if needed
sudo apt install ghostscript

# Convert to PNG
gs -sDEVICE=png16m -sOutputFile=structure.png -r300 -dNOPAUSE -dBATCH rna_ss.ps
gs -sDEVICE=png16m -sOutputFile=dotplot.png -r300 -dNOPAUSE -dBATCH rna_dp.ps

# Convert to PDF
gs -sDEVICE=pdfwrite -sOutputFile=structure.pdf -dNOPAUSE -dBATCH rna_ss.ps
gs -sDEVICE=pdfwrite -sOutputFile=dotplot.pdf -dNOPAUSE -dBATCH rna_dp.ps
```

## Useful RNAfold Options

| Option | Description |
|--------|-------------|
| `-p` | Calculate partition function + base pair probabilities |
| `-T <temp>` | Temperature in Celsius (default: 37) |
| `--MEA` | Calculate Maximum Expected Accuracy structure |
| `-d2` | Allow dangles on both sides of helices |
| `--noPS` | Don't generate structure PostScript |
| `--noDP` | Don't generate dot plot PostScript |
| `-C` | Read structure constraints from input |

## Using with Your Generated Sequences

```bash
# Extract sequences and process them
for i in {1..5}; do
    SEQUENCE=$(grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand${i}.txt | cut -d' ' -f2)
    echo ">candidate_${i}" > candidate_${i}.fa
    echo "$SEQUENCE" >> candidate_${i}.fa
    RNAfold -p < candidate_${i}.fa
done
```

## Batch Processing Multiple Candidates
```bash
# Create one file with all candidates
cat > all_candidates.fa << EOF
$(for i in {1..5}; do
    echo ">candidate_$i"
    grep "Sequence:" designed_sequences/z_tile_tetramer_no_filter_cand${i}.txt | cut -d' ' -f2
done)
EOF

# Process all at once
RNAfold -p < all_candidates.fa
```

This gives you the exact same results as the web interface, with PostScript files you can view or convert to PNG/PDF!