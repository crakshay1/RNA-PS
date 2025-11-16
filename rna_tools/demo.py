# =========================
# Demo
# =========================

import create_rna_data as cr
import rna_visualizer as rv
import process_rna_data as prd

if __name__ == "__main__":
    scaffold_name = cr.eg.choicebox("Pick one of the sccafolds", "Scaffold", cr.scaffold.keys())
    cands = cr.generate_candidates_for_scaffold(scaffold_name, n_candidates=cr.eg.integerbox("Number of candidates"))

    cr.save_candidates(scaffold_name, cands)

    # Designing
    count = 0
    for candidate in cands:
        count += 1
        seq = candidate["sequence"]
        ss = candidate["predicted_ss"]  # or best["target_ss"]
        # 1) Write inputs for RNAComposer
        rna_composer_input = cr.write_rnacomposer_input(seq, ss)
        # 2) Creates the .pdb file
        cr.create_pdb_from_RNAComposer(rna_composer_input, count)
    
    # Checking RNA stability
    prd.os.makedirs("MFE_test", exist_ok=True)
    prd.os.path.join("MFE_test", f'output{count}')
    dict_MFE = {}
    for file in rv.os.listdir("pdb_files"):
        prd.process_structure_file(f"pdb_files/{file}", f"MFE_test/output_{rv.os.path.basename(file)}")
        dict_MFE[file] = 0
    with open("analysis/energy.txt") as f:
        lines = f.read().splitlines()
    prd.os.remove("analysis/energy.txt")
    for i, key in enumerate(dict_MFE):
        dict_MFE[key] = float(lines[i])
    dict_MFE_sorted = dict(sorted(dict_MFE.items(), key=lambda item: item[1]))
    top1 = list(dict_MFE_sorted.keys())[0]

    # Displaying RNA structure in an interactive window
    rv.represent(f"pdb_files/{top1}")