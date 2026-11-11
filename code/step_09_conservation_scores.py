#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 9: Conservation Score Per Position")
    print("=" * 65)

    msa_path = os.path.join(DATA_RAW, "p53_msa.fasta")
    if not os.path.exists(msa_path):
        print("[ERROR] Run step_08 first!")
        sys.exit(1)

    # Parse MSA
    sequences = {}
    current = None
    with open(msa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:]
                sequences[current] = ""
            elif current:
                sequences[current] += line

    print("[Step 9] Sequences in MSA: " + str(len(sequences)))

    # Find human sequence
    human_key = None
    for k in sequences:
        if "Human" in k or "human" in k or "Homo" in k:
            human_key = k
            break
    if not human_key:
        human_key = list(sequences.keys())[0]
    print("[Step 9] Human sequence key: " + human_key)

    human_aligned = sequences[human_key]
    alignment_len = len(human_aligned)
    print("[Step 9] Alignment length: " + str(alignment_len) + " columns")

    # Load WT sequence for verification
    with open(os.path.join(DATA_RAW, "p53_sequence.txt")) as f:
        wt_seq = f.read().strip()

    # Map alignment columns to protein positions
    # Human sequence has gaps in alignment; skip gap columns
    col_to_pos = {}
    pos = 0
    for col in range(alignment_len):
        if human_aligned[col] != "-":
            pos += 1
            col_to_pos[col] = pos

    print("[Step 9] Human residues mapped: " + str(pos) + " positions")

    # Calculate conservation for each position
    all_seqs = list(sequences.values())
    num_seqs = len(all_seqs)

    conservation = {}
    for col, protein_pos in col_to_pos.items():
        human_aa = human_aligned[col]
        match_count = 0
        total_count = 0
        for seq in all_seqs:
            if col < len(seq):
                aa = seq[col]
                if aa != "-":
                    total_count += 1
                    if aa == human_aa:
                        match_count += 1

        score = match_count / total_count if total_count > 0 else 0
        conservation[protein_pos] = round(score, 4)

    print("[Step 9] Conservation scores computed for " + str(len(conservation)) + " positions")

    # Save
    out_path = os.path.join(DATA_INTER, "conservation_scores.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["position", "wt_aa", "conservation_score"])
        for p in range(1, SEQUENCE_LEN + 1):
            aa = wt_seq[p - 1] if p <= len(wt_seq) else "?"
            score = conservation.get(p, 0)
            w.writerow([p, aa, score])

    scores = [conservation.get(p, 0) for p in range(1, SEQUENCE_LEN + 1)]
    avg = round(sum(scores) / len(scores), 4)

    # Domain averages
    print("[Step 9] Conservation by domain:")
    for code, (start, end, name) in DOMAINS.items():
        domain_scores = [conservation.get(p, 0) for p in range(start, end + 1)]
        davg = round(sum(domain_scores) / len(domain_scores), 4) if domain_scores else 0
        print("    " + code + " (" + str(start) + "-" + str(end) + "): " + str(davg))

    # Zinc residue conservation
    print("[Step 9] Zinc residue conservation:")
    for pos in sorted(ZINC_RESIDUES.keys()):
        score = conservation.get(pos, 0)
        print("    " + ZINC_RESIDUES[pos] + str(pos) + ": " + str(score))

    print("")
    print("=" * 65)
    print("  STEP 9 COMPLETE")
    print("=" * 65)
    print("  Positions scored:  " + str(len(conservation)) + " / " + str(SEQUENCE_LEN))
    print("  Average conserv:   " + str(avg))
    print("  Most conserved:    " + str(max(scores)))
    print("  Least conserved:   " + str(min(scores)))
    print("  Output: " + out_path + " (" + str(os.path.getsize(out_path)) + " bytes)")
    print("")
    print("  Next: python code/step_10_physicochemical.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
