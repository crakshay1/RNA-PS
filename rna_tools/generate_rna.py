
import os
import random
import RNA  # ViennaRNA Python bindings

# =========================
# 1. Scaffolds (yours)
# =========================

scaffold = {
    "z_tile_tetramer": "((((....))))..((((....))))..((((....))))..((((....))))",
    "tetrahedron_wireframe": "((((....))))..((((....))))..((((....))))..((((....))))",
    "triangular_prism": (
        "((((....))))..((((....))))..((((....))))"
        "........"
        "((((....))))..((((....))))..((((....))))"
    ),
    "rool_repeat_unit": "(((....)))........(((....)))........(((....)))........(((....)))",
}

NUCS = ["A", "U", "G", "C"]

# =========================
# 2. Loop detection
# =========================

def find_hairpin_loops(dot_bracket: str):
    """
    Return a list of hairpin loops as (start, end) indices.
    Hairpin = run of '.' that is immediately preceded by '(' and followed by ')'.
    """
    loops = []
    i = 0
    n = len(dot_bracket)
    while i < n:
        if dot_bracket[i] == '.':
            j = i
            while j < n and dot_bracket[j] == '.':
                j += 1
            # i..j-1 is a run of '.'
            if i > 0 and j < n and dot_bracket[i-1] == '(' and dot_bracket[j] == ')':
                loops.append((i, j-1))
            i = j
        else:
            i += 1
    return loops

# =========================
# 3. Motif config sampling
# =========================

def sample_motif_configuration(scaffold_name: str, loops: list, rng=None):
    """
    Choose which loops get which motif type.
    Returns: dict keyed by loop index (0,1,2,...) with entries like:
      {"type": "GNRA"} or {"type": "kissing", "pair_with": 2}
    """
    if rng is None:
        rng = random.Random()

    n_loops = len(loops)
    config = {i: {"type": None} for i in range(n_loops)}

    if scaffold_name == "z_tile_tetramer":
        # 4 loops: choose a kissing pairing pattern among a few options
        if n_loops != 4:
            return config  # fallback

        patterns = [
            [(0, 1), (2, 3)],
            [(0, 2), (1, 3)],
            [(0, 3), (1, 2)],
        ]
        pairs = rng.choice(patterns)
        for a, b in pairs:
            config[a] = {"type": "kissing", "pair_with": b}
            config[b] = {"type": "kissing", "pair_with": a}

    elif scaffold_name == "tetrahedron_wireframe":
        # Opposite loops kissing: (0<->2, 1<->3)
        if n_loops >= 4:
            config[0] = {"type": "kissing", "pair_with": 2}
            config[2] = {"type": "kissing", "pair_with": 0}
            config[1] = {"type": "kissing", "pair_with": 3}
            config[3] = {"type": "kissing", "pair_with": 1}

    elif scaffold_name == "triangular_prism":
        # Top (0,1,2) vs bottom (3,4,5) loops, pair them if present
        if n_loops >= 6:
            for top, bottom in zip(range(0, 3), range(3, 6)):
                config[top] = {"type": "kissing", "pair_with": bottom}
                config[bottom] = {"type": "kissing", "pair_with": top}

    # For any loop not assigned kissing, assign a stabilizing tetraloop motif
    for i in range(n_loops):
        if config[i]["type"] is None:
            # randomly choose GNRA or UUCG
            config[i] = {"type": random.choice(["GNRA", "UUCG"])}

    return config

# =========================
# 4. Motifs -> constraints + annotation
# =========================

def random_base(rng=None):
    if rng is None:
        rng = random.Random()
    return rng.choice(NUCS)

def complement(base: str) -> str:
    return {"A": "U", "U": "A", "G": "C", "C": "G"}.get(base, "A")

def motifs_to_constraints(dot_bracket: str, loops: list, motif_config: dict, rng=None):
    """
    Convert motif configuration into:
      - constraints: {position: [allowed_bases]}
      - annotation: list of same length as dot_bracket, with motif labels or None
    """
    if rng is None:
        rng = random.Random()

    n = len(dot_bracket)
    constraints = {}
    annotation = [None] * n

    # Helper to set allowed bases
    def set_allowed(pos, bases):
        constraints[pos] = list(bases)

    # First handle non-kissing motifs (GNRA, UUCG)
    for loop_idx, (start, end) in enumerate(loops):
        motif = motif_config[loop_idx]
        mtype = motif["type"]

        if mtype == "GNRA":
            length = end - start + 1
            if length < 4:
                continue
            # Use the first 4 positions of the loop
            i0 = start
            # G
            set_allowed(i0, ["G"])
            # N
            set_allowed(i0 + 1, NUCS)
            # R (A/G)
            set_allowed(i0 + 2, ["A", "G"])
            # A
            set_allowed(i0 + 3, ["A"])
            for k in range(4):
                annotation[i0 + k] = "GNRA"

        elif mtype == "UUCG":
            length = end - start + 1
            if length < 4:
                continue
            i0 = start
            seq = ["U", "U", "C", "G"]
            for k, b in enumerate(seq):
                set_allowed(i0 + k, [b])
                annotation[i0 + k] = "UUCG"

    # Then handle kissing-loop pairs
    done_pairs = set()
    for loop_idx, (start, end) in enumerate(loops):
        motif = motif_config[loop_idx]
        if motif["type"] != "kissing":
            continue
        pair_with = motif.get("pair_with")
        if pair_with is None:
            continue
        # Avoid double-processing same pair
        key = tuple(sorted((loop_idx, pair_with)))
        if key in done_pairs:
            continue
        done_pairs.add(key)

        # This loop and partner loop
        (s1, e1) = loops[loop_idx]
        (s2, e2) = loops[pair_with]

        len1 = e1 - s1 + 1
        len2 = e2 - s2 + 1
        L = min(len1, len2)  # use min length

        # Generate random sequence for loop1, enforce reverse complement on loop2
        seq1 = [random_base(rng) for _ in range(L)]
        seq2 = [complement(b) for b in reversed(seq1)]

        for k in range(L):
            pos1 = s1 + k
            pos2 = s2 + k  # aligned positions; other alignments also possible

            set_allowed(pos1, [seq1[k]])
            set_allowed(pos2, [seq2[k]])

            annotation[pos1] = f"KL_{loop_idx}_{pair_with}"
            annotation[pos2] = f"KL_{pair_with}_{loop_idx}"

    return constraints, annotation


# =========================
# 5. 
# =========================

def build_initial_sequence(length, constraints, rng=None):
    """
    Build a random initial RNA sequence of given length that respects
    the base constraints at certain positions.
    constraints: dict {pos: [allowed_bases]}
    """
    if rng is None:
        rng = random.Random()
    seq = []
    for i in range(length):
        if i in constraints:
            allowed = constraints[i]
            seq.append(rng.choice(allowed))
        else:
            seq.append(rng.choice(NUCS))
    return "".join(seq)


# =========================
# 6. 
# =========================

def inverse_fold_with_constraints(target_ss, constraints, n_tries=10, rng=None):
    """
    Try multiple times to design a sequence whose MFE structure matches target_ss
    as closely as possible, given base constraints.
    Returns: best (sequence, predicted_ss, mfe, bp_distance)
    """
    if rng is None:
        rng = random.Random()

    best = None
    best_dist = None

    for _ in range(n_tries):
        # 1) build an initial sequence that respects constraints
        init_seq = build_initial_sequence(len(target_ss), constraints, rng)

        # 2) run inverse folding from ViennaRNA
        # NOTE: depending on your ViennaRNA version, inverse_fold signature may differ.
        designed_seq, _ = RNA.inverse_fold(init_seq, target_ss)

        # 3) refold the designed sequence to see what it actually does
        pred_ss, mfe = RNA.fold(designed_seq)
        dist = RNA.bp_distance(target_ss, pred_ss)

        if best is None or dist < best_dist:
            best = (designed_seq, pred_ss, mfe, dist)
            best_dist = dist

        # early stop if we get a perfect match
        if dist == 0:
            break

    return best  # (sequence, predicted_ss, mfe, bp_distance)


# =========================
# 7. 
# =========================

def generate_candidates_for_scaffold(
    scaffold_name: str,
    n_candidates: int = 5,
    rng_seed: int = 42,
):
    """
    Full pipeline for one scaffold:
      - get dot-bracket
      - find loops
      - sample motif configuration
      - convert motifs -> constraints
      - run inverse folding multiple times
    Returns a list of dicts with sequence, structure, mfe, etc.
    """
    rng = random.Random(rng_seed)
    db = scaffold[scaffold_name]

    # 1) find hairpin loops in the scaffold
    loops = find_hairpin_loops(db)

    # 2) sample motif configuration for this scaffold
    motif_cfg = sample_motif_configuration(scaffold_name, loops, rng)

    # 3) convert motifs -> constraints + annotation
    constraints, annotation = motifs_to_constraints(db, loops, motif_cfg, rng)

    candidates = []

    # oversample a bit so we can pick the best n_candidates
    n_trials = max(n_candidates * 3, n_candidates)

    for t in range(n_trials):
        result = inverse_fold_with_constraints(db, constraints, n_tries=5, rng=rng)
        if result is None:
            continue
        seq, pred_ss, mfe, dist = result
        candidates.append({
            "sequence": seq,
            "predicted_ss": pred_ss,
            "target_ss": db,
            "mfe": mfe,
            "bp_distance": dist,
            "motifs": motif_cfg,
            "annotation": annotation,
        })

    # sort by base-pair distance, then by MFE (more negative = better)
    candidates.sort(key=lambda c: (c["bp_distance"], c["mfe"]))

    return candidates[:n_candidates]


# =========================
# 8. Save generated candidates
# =========================

def save_candidates(scaffold_name, candidates, out_dir="designed_sequences"):
    os.makedirs(out_dir, exist_ok=True)
    for i, c in enumerate(candidates):
        base = f"{scaffold_name}_cand{i+1}"
        txt_path = os.path.join(out_dir, base + ".txt")
        with open(txt_path, "w") as f:
            f.write(f"Name: {base}\n")
            f.write(f"Sequence: {c['sequence']}\n")
            f.write(f"Target_SS: {c['target_ss']}\n")
            f.write(f"Predicted_SS: {c['predicted_ss']}\n")
            f.write(f"MFE: {c['mfe']}\n")
            f.write(f"BP_distance: {c['bp_distance']}\n")
