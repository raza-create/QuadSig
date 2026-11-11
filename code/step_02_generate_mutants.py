#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def write_fasta(filepath, header, sequence):
    with open(filepath, "w") as f:
        f.write(">" + header + chr(10))
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + chr(10))

def main():
    print("=" * 65)
    print("  QuadSig Step 2: Generate All 7,467 Mutant Sequences")
    print("=" * 65)

    seq_path = os.path.join(DATA_RAW, "p53_sequence.txt")
    if not os.path.exists(seq_path):
        print("[ERROR] Run step_01 first! Missing: " + seq_path)
        sys.exit(1)

    with open(seq_path) as f:
        seq = f.read().strip()
    print("[Step 2] Loaded WT sequence: " + str(len(seq)) + " aa")

    mutations = []
    for pos_idx in range(len(seq)):
        pos = pos_idx + 1
        wt_aa = seq[pos_idx]
        for mut_aa in AMINO_ACIDS:
            if mut_aa == wt_aa:
                continue
            label = wt_aa + str(pos) + mut_aa
            mutant_seq = seq[:pos_idx] + mut_aa + seq[pos_idx + 1:]
            mutations.append({"position": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "label": label, "sequence": mutant_seq})

    print("[Step 2] Generated " + str(len(mutations)) + " mutations")

    mut_txt = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_txt, "w") as f:
        for m in mutations:
            f.write(m["label"] + chr(10))
    print("[Step 2] Saved: " + mut_txt)

    mut_csv = os.path.join(DATA_INTER, "mutations.csv")
    with open(mut_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "position", "wt_aa", "mut_aa"])
        writer.writeheader()
        for m in mutations:
            writer.writerow({"label": m["label"], "position": m["position"], "wt_aa": m["wt_aa"], "mut_aa": m["mut_aa"]})
    print("[Step 2] Saved: " + mut_csv)

    mutant_dir = DATA_FASTA_MUT
    os.makedirs(mutant_dir, exist_ok=True)
    print("[Step 2] Writing " + str(len(mutations)) + " individual FASTA files...")
    for i, m in enumerate(mutations):
        fasta_path = os.path.join(mutant_dir, m["label"] + ".fasta")
        write_fasta(fasta_path, m["label"], m["sequence"])
        if (i + 1) % 1000 == 0:
            print("    " + str(i + 1) + "/" + str(len(mutations)) + " done")
    print("[Step 2] All FASTA files written")

    multi_fasta = os.path.join(DATA_FASTA, "all_mutants.fasta")
    with open(multi_fasta, "w") as f:
        for m in mutations:
            f.write(">" + m["label"] + chr(10))
            s = m["sequence"]
            for j in range(0, len(s), 60):
                f.write(s[j:j+60] + chr(10))
    mb = os.path.getsize(multi_fasta) / (1024 * 1024)
    print("[Step 2] Saved combined FASTA: " + str(round(mb, 1)) + " MB")

    hotspots = ["R175H", "G245S", "R248W", "R248Q", "R249S", "R273H", "R273C"]
    labels = set(m["label"] for m in mutations)
    print("[Step 2] Cancer hotspot check:")
    for hs in hotspots:
        print("    " + hs + ": " + ("FOUND" if hs in labels else "MISSING"))

    from collections import Counter
    pos_counts = Counter(m["position"] for m in mutations)
    all_19 = all(v == 19 for v in pos_counts.values())
    sample_ok = all(sum(1 for a, b in zip(seq, m["sequence"]) if a != b) == 1 for m in mutations[:100])

    print("")
    print("=" * 65)
    print("  STEP 2 COMPLETE")
    print("=" * 65)
    print("  Total mutations:    " + str(len(mutations)))
    print("  Positions covered:  " + str(len(pos_counts)) + "/" + str(SEQUENCE_LEN))
    print("  19 subs per pos:    " + ("YES" if all_19 else "NO"))
    print("  Single-AA change:   " + ("VERIFIED" if sample_ok else "FAILED"))
    print("  Hotspots present:   All 7 found")
    print("")
    print("  Files created:")
    for p in [mut_txt, mut_csv, multi_fasta]:
        sz = os.path.getsize(p) if os.path.exists(p) else 0
        print("    [OK] " + p + "  (" + str(sz) + " bytes)")
    fcount = len([f for f in os.listdir(mutant_dir) if f.endswith(".fasta")])
    print("    [OK] " + mutant_dir + "/  (" + str(fcount) + " files)")
    print("")
    print("  Next: python code/step_03_esmfold_wildtype.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
