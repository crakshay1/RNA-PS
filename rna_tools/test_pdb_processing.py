#!/usr/bin/env python3

import os
import process_rna_data as prd

def test_pdb_processing():
    """Test PDB file processing capabilities."""

    print("=== Testing PDB File Processing ===")

    # Create a simple test PDB file (minimal RNA structure)
    test_pdb_content = """HEADER    RNA                             16-NOV-25   TEST
ATOM      1  P     A A   1      -0.123   1.456   0.987  1.00 20.00           P
ATOM      2  O1P   A A   1      -1.234   2.567   1.098  1.00 20.00           O
ATOM      3  O2P   A A   1       0.890   1.345   2.109  1.00 20.00           O
ATOM      4  O5'   A A   1       0.456   0.123  -0.234  1.00 20.00           O
ATOM      5  C5'   A A   1       1.567   0.234  -1.345  1.00 20.00           C
ATOM      6  C4'   A A   1       2.678  -0.890  -1.456  1.00 20.00           C
ATOM      7  O4'   A A   1       3.789  -0.789  -0.567  1.00 20.00           O
ATOM      8  C3'   A A   1       3.012  -1.123  -2.890  1.00 20.00           C
ATOM      9  O3'   A A   1       2.234  -2.234  -3.456  1.00 20.00           O
ATOM     10  C2'   A A   1       4.456  -1.567  -2.678  1.00 20.00           C
ATOM     11  O2'   A A   1       5.123  -1.890  -3.890  1.00 20.00           O
ATOM     12  C1'   A A   1       4.890  -0.456  -1.789  1.00 20.00           C
ATOM     13  N9    A A   1       5.567   0.567  -2.567  1.00 20.00           N
ATOM     14  C8    A A   1       5.890   1.890  -2.234  1.00 20.00           C
ATOM     15  N7    A A   1       6.567   2.567  -3.012  1.00 20.00           N
ATOM     16  C5    A A   1       6.789   1.789  -4.123  1.00 20.00           C
ATOM     17  C6    A A   1       7.456   1.890  -5.345  1.00 20.00           C
ATOM     18  N6    A A   1       8.123   2.890  -5.678  1.00 20.00           N
ATOM     19  N1    A A   1       7.345   0.789  -6.123  1.00 20.00           N
ATOM     20  C2    A A   1       6.567  -0.234  -5.789  1.00 20.00           C
ATOM     21  N3    A A   1       5.890  -0.345  -4.678  1.00 20.00           N
ATOM     22  C4    A A   1       6.012   0.567  -3.890  1.00 20.00           C
ATOM     23  P     U A   2       1.567  -3.456  -2.890  1.00 20.00           P
ATOM     24  O1P   U A   2       0.456  -4.567  -3.567  1.00 20.00           O
ATOM     25  O2P   U A   2       2.678  -4.123  -2.234  1.00 20.00           O
ATOM     26  O5'   U A   2       0.890  -2.890  -1.567  1.00 20.00           O
ATOM     27  C5'   U A   2      -0.234  -3.567  -0.890  1.00 20.00           C
ATOM     28  C4'   U A   2      -1.345  -2.890  -0.123  1.00 20.00           C
ATOM     29  O4'   U A   2      -2.456  -2.567   0.789  1.00 20.00           O
ATOM     30  C3'   U A   2      -0.890  -1.567   0.567  1.00 20.00           C
ATOM     31  O3'   U A   2      -1.567  -0.456   1.234  1.00 20.00           O
ATOM     32  C2'   U A   2      -2.012  -1.890   1.567  1.00 20.00           C
ATOM     33  O2'   U A   2      -3.234  -1.234   1.890  1.00 20.00           O
ATOM     34  C1'   U A   2      -2.567  -3.234   0.890  1.00 20.00           C
ATOM     35  N1    U A   2      -3.890  -3.567   1.567  1.00 20.00           N
ATOM     36  C2    U A   2      -4.567  -4.789   1.234  1.00 20.00           C
ATOM     37  O2    U A   2      -4.234  -5.456   0.234  1.00 20.00           O
ATOM     38  N3    U A   2      -5.789  -5.012   2.012  1.00 20.00           N
ATOM     39  C4    U A   2      -6.234  -4.345   3.123  1.00 20.00           C
ATOM     40  O4    U A   2      -7.345  -4.678   3.789  1.00 20.00           O
ATOM     41  C5    U A   2      -5.456  -3.123   3.456  1.00 20.00           C
ATOM     42  C6    U A   2      -4.234  -2.890   2.678  1.00 20.00           C
TER      43        U A   2
END                                                                             """

    # Write test PDB file
    test_pdb_file = "test_rna.pdb"
    with open(test_pdb_file, "w") as f:
        f.write(test_pdb_content)

    print(f"Created test PDB file: {test_pdb_file}")

    # Test processing
    try:
        print("\\nTesting PDB processing...")
        prd.process_structure_file(test_pdb_file, "pdb_output")

        print("\\nTesting specific PDB function...")
        prd.process_pdb_file(test_pdb_file, "pdb_output_2")

        print("\\nTesting with non-existent file...")
        prd.process_structure_file("non_existent.pdb", "test_output")

        print("\\nTesting with unsupported format...")
        with open("test.txt", "w") as f:
            f.write("This is not a structure file")
        prd.process_structure_file("test.txt", "test_output")

    except Exception as e:
        print(f"Error during testing: {e}")

    finally:
        # Cleanup
        for file in [test_pdb_file, "test.txt"]:
            if os.path.exists(file):
                os.remove(file)

    print("\\n=== PDB Processing Test Complete ===")

if __name__ == "__main__":
    test_pdb_processing()