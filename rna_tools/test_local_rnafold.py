#!/usr/bin/env python3

import os
import re
import subprocess

def extract_sequence_from_file(file_path):
    """Extract RNA sequence from candidate file."""
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('Sequence:'):
                return line.split(':', 1)[1].strip()
    return None

def test_local_rnafold():
    """Test local RNAfold with existing sequence files."""

    # Find the best candidate files (usually have lowest MFE)
    sequence_dir = "designed_sequences"
    test_files = [
        "z_tile_tetramer_no_filter_cand1.txt",
        "z_tile_tetramer_percentile_cand1.txt",
        "z_tile_tetramer_threshold_cand1.txt"
    ]

    print("=== Testing Local RNAfold with Your Generated Sequences ===\n")

    for i, test_file in enumerate(test_files, 1):
        file_path = os.path.join(sequence_dir, test_file)

        if not os.path.exists(file_path):
            print(f"Skipping {test_file} - file not found")
            continue

        sequence = extract_sequence_from_file(file_path)
        if not sequence:
            print(f"Skipping {test_file} - no sequence found")
            continue

        print(f"{i}. Testing with {test_file}")
        print(f"   Sequence: {sequence}")
        print(f"   Length: {len(sequence)} nucleotides")

        # Test 1: Basic folding (no plot)
        print("   Running basic fold...")
        cmd = ["python", "local_rnafold.py", "-s", sequence, "--no-plot"]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Extract MFE from output
            mfe_match = re.search(r'MFE Energy: ([-\d.]+)', result.stdout)
            if mfe_match:
                mfe = float(mfe_match.group(1))
                print(f"   ✓ MFE: {mfe:.2f} kcal/mol")

        except subprocess.CalledProcessError as e:
            print(f"   ✗ Error: {e}")
            continue

        # Test 2: Generate plot
        output_png = f"rnafold_test_{i}.png"
        print(f"   Generating plot: {output_png}")
        cmd = ["python", "local_rnafold.py", "-s", sequence, "-o", output_png]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"   ✓ Plot saved to {output_png}")
        except subprocess.CalledProcessError:
            print(f"   ✗ Plot generation failed")

        print()

    # Test 3: Test with different temperatures
    if test_files and os.path.exists(os.path.join(sequence_dir, test_files[0])):
        print("Testing temperature effects...")
        sequence = extract_sequence_from_file(os.path.join(sequence_dir, test_files[0]))

        for temp in [25, 37, 60]:
            print(f"   Temperature {temp}°C...")
            cmd = ["python", "local_rnafold.py", "-s", sequence, "-t", str(temp), "--no-plot"]
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                mfe_match = re.search(r'MFE Energy: ([-\d.]+)', result.stdout)
                if mfe_match:
                    mfe = float(mfe_match.group(1))
                    print(f"     MFE: {mfe:.2f} kcal/mol")
            except subprocess.CalledProcessError:
                print(f"     Error at {temp}°C")

        print()

    print("=== Test Commands You Can Run ===")
    print("1. Basic folding:")
    print("   source .venv/bin/activate && python local_rnafold.py -s CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC")
    print()
    print("2. With visualization:")
    print("   source .venv/bin/activate && python local_rnafold.py -s CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC -o my_result.png")
    print()
    print("3. Different temperature:")
    print("   source .venv/bin/activate && python local_rnafold.py -s CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC -t 25")
    print()
    print("4. From file (create sequence.txt first):")
    print("   echo 'CCUCAAGUGAGGCACCCCUUAAGGGGACGAUAUUAAUGUCACGCACACUUGUGC' > sequence.txt")
    print("   source .venv/bin/activate && python local_rnafold.py -f sequence.txt")

if __name__ == "__main__":
    test_local_rnafold()