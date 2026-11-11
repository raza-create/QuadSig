#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 12: Zinc-Coordination Residue Flagging")
    print("=" * 65)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 12] Total mutations: " + str(len(all_mutations)))
    print("[Step 12] Zinc residues: " + str(ZINC_RESIDUES))

    # Load structural scores for validation check
    struct_path = os.path.join(DATA_INTER, "structural_scores.csv")
    struct_data = {}
    if os.path.exists(struct_path):
        with open(struct_path) as f:
            for row in csv.DictReader(f):
                struct_data[row["mutation"]] = row

    results = []
    zinc_count = 0
    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        is_zinc = pos in ZINC_RESIDUES

        if is_zinc:
            zinc_count += 1
            label = ZINC_RESIDUES[pos] + str(pos) + " zinc coordination"
        else:
            label = ""

        results.append({
            "mutation": mutation,
            "position": pos,
            "wt_aa": mutation[0],
            "mut_aa": mutation[-1],
            "is_zinc_coord": is_zinc,
            "zinc_position_label": label,
            "expected_score_threshold": 0.90 if is_zinc else "",
        })

    out_path = os.path.join(DATA_INTER, "zinc_flags.csv")
    fields = ["mutation","position","wt_aa","mut_aa","is_zinc_coord","zinc_position_label","expected_score_threshold"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    print("[Step 12] Zinc mutations flagged: " + str(zinc_count))
    print("[Step 12] Non-zinc mutations: " + str(len(results) - zinc_count))

    # Validation: check structural scores for zinc mutations
    print("[Step 12] Zinc mutation structural impact:")
    for pos in sorted(ZINC_RESIDUES.keys()):
        aa = ZINC_RESIDUES[pos]
        zinc_muts = [r for r in results if r["position"] == pos and r["is_zinc_coord"]]
        for zm in zinc_muts:
            sd = struct_data.get(zm["mutation"], {})
            rmsd = sd.get("rmsd", "?")
            plddt_change = sd.get("plddt_change", "?")
            print("    " + zm["mutation"] + ": RMSD=" + str(rmsd) + " pLDDT_change=" + str(plddt_change))

    print("")
    print("=" * 65)
    print("  STEP 12 COMPLETE")
    print("=" * 65)
    print("  Zinc mutations:    " + str(zinc_count) + " (4 positions x 19 subs = 76)")
    print("  Validation:        All zinc subs should score > 0.90 in final QuadSig")
    print("  Output: " + out_path)
    print("  Next: python code/step_13_gof_detection.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
