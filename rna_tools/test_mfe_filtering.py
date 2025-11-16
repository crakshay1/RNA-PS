#!/usr/bin/env python3

import generate_rna as gr

def test_mfe_filtering():
    """Test MFE filtering functionality with different parameters."""

    print("=== Testing MFE Filtering ===")
    scaffold_name = "z_tile_tetramer"

    print(f"\n1. Generating candidates without MFE filtering...")
    candidates_no_filter = gr.generate_candidates_for_scaffold(
        scaffold_name,
        n_candidates=15, 
        rng_seed=42
    )

    if candidates_no_filter:
        stats = gr.get_mfe_statistics(candidates_no_filter)
        print(f"Generated {len(candidates_no_filter)} candidates")
        print(f"MFE range: {stats['min_mfe']:.2f} to {stats['max_mfe']:.2f}")
        print(f"Mean MFE: {stats['mean_mfe']:.2f}")

        print(f"\n2. Testing MFE threshold filtering...")
        threshold = stats['mean_mfe']  
        candidates_threshold = gr.generate_candidates_for_scaffold(
            scaffold_name,
            n_candidates=15,
            rng_seed=42,
            mfe_threshold=threshold
        )

        print(f"\n3. Testing MFE percentile filtering (bottom 30%)...")
        candidates_percentile = gr.generate_candidates_for_scaffold(
            scaffold_name,
            n_candidates=15,
            rng_seed=42,
            mfe_percentile=30 
        )

        print(f"\n4. Testing combined MFE filtering...")
        candidates_combined = gr.generate_candidates_for_scaffold(
            scaffold_name,
            n_candidates=15,
            rng_seed=42,
            mfe_threshold=stats['max_mfe'],  
            mfe_percentile=50 
        )

        print(f"\n5. Saving results...")
        gr.save_candidates(f"{scaffold_name}_no_filter", candidates_no_filter)
        gr.save_candidates(f"{scaffold_name}_threshold", candidates_threshold)
        gr.save_candidates(f"{scaffold_name}_percentile", candidates_percentile)
        gr.save_candidates(f"{scaffold_name}_combined", candidates_combined)

        print("\nResults saved to designed_sequences/ directory")
        print("Check the summary files to compare MFE statistics!")

    else:
        print("No candidates generated")

if __name__ == "__main__":
    test_mfe_filtering()