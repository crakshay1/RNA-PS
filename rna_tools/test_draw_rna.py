
import os
from draw_rna import draw
from generate_rna import scaffold

def test_draw_rna():
    sequence = "GGGGAAAACCCC"
    structure = "((((....))))"
    output_png_path = "test_draw_rna.png"

    try:
        draw.draw_rna(sequence, structure, output_png_path)
        print(f"Test visualization saved to {output_png_path}")
    except Exception as e:
        print(f"Error generating test visualization: {e}")

def test_draw_scaffold():
    

if __name__ == "__main__":
    test_draw_rna()
