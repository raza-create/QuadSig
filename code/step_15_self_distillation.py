#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 15: Self-Distillation Round 2 Label Cleaning")
    print("=" * 65)

    scores_path = os.path.join(DATA_INTER, "quadsig_scores.csv")
    with open(scores_path) as f:
        scores = {r["mutation"]: r for r in csv.DictReader(f)}
    print("[Step 15] Round 1 scores loaded: " + str(len(scores)))

    pop = {}
    with open(os.path.join(DATA_INTER, "population_scores.csv")) as f:
        for r in csv.DictReader(f):
            pop[r["mutation"]] = r

    cons = {}
    with open(os.path.join(DATA_INTER, "conservation_scores.csv")) as f:
        for r in csv.DictReader(f):
            cons[int(r["position"])] = float(r["conservation_score"])

    struct = {}
    with open(os.path.join(DATA_INTER, "structural_scores.csv")) as f:
        for r in csv.DictReader(f):
            struct[r["mutation"]] = r

    with open(os.path.join(DATA_INTER, "mutations.txt")) as f:
        all_mutations = [line.strip() for line in f if line.strip()]

    removed = []
    cleaned = []

    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        s = scores.get(mutation, {})
        r1_score = float(s.get("quadsig_score", 0.5))
        p = pop.get(mutation, {})
        pop_score = float(p.get("population_score", 1.0))
        conservation = cons.get(pos, 0.5)
        sd = struct.get(mutation, {})
        rmsd = float(sd.get("rmsd", 0))

        is_noisy = (pop_score >= 0.95 and r1_score < 0.20)

        if is_noisy:
            cleaned.append({
                "mutation": mutation,
                "round1_score": r1_score,
                "original_label": s.get("verdict", "uncertain"),
                "cleaned_label": "UNCERTAIN",
                "reason": "Low structural/evolutionary signal despite population rarity",
                "conservation": round(conservation, 4),
                "rmsd": rmsd,
            })
            removed.append(mutation)
        else:
            cleaned.append({
                "mutation": mutation,
                "round1_score": r1_score,
                "original_label": s.get("verdict", "uncertain"),
                "cleaned_label": s.get("verdict", "uncertain"),
                "reason": "",
                "conservation": round(conservation, 4),
                "rmsd": rmsd,
            })

    out_path = os.path.join(DATA_INTER, "cleaned_labels.csv")
    fields = ["mutation","round1_score","original_label","cleaned_label","reason","conservation","rmsd"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(cleaned)

    domain_removed = {}
    for mut in removed:
        pos = int(mut[1:-1])
        found = "LINKER"
        for code, (start, end, name) in DOMAINS.items():
            if start <= pos <= end:
                found = code
                break
        domain_removed[found] = domain_removed.get(found, 0) + 1

    print("[Step 15] Removed from pathogenic group: " + str(len(removed)))
    print("[Step 15] Relabeled as UNCERTAIN (not BENIGN)")
    print("[Step 15] Removal by domain:")
    for d, count in sorted(domain_removed.items(), key=lambda x: -x[1]):
        print("    " + d + ": " + str(count))

    hotspots = ["R175H","G245S","R248W","R248Q","R249S","R273H","R273C"]
    print("[Step 15] Cancer hotspot check:")
    for hs in hotspots:
        status = "REMOVED" if hs in removed else "KEPT"
        print("    " + hs + ": " + status)

    for c in cleaned:
        if c["mutation"] == "P72R":
            print("[Step 15] P72R: label=" + c["cleaned_label"] + " score=" + str(c["round1_score"]))

    final_harmful = sum(1 for c in cleaned if c["cleaned_label"] == "harmful")
    final_uncertain = sum(1 for c in cleaned if c["cleaned_label"] in ["uncertain", "UNCERTAIN"])
    final_benign = sum(1 for c in cleaned if c["cleaned_label"] == "benign")

    print("")
    print("=" * 65)
    print("  STEP 15 COMPLETE")
    print("=" * 65)
    print("  Total mutations:     " + str(len(cleaned)))
    print("  Removed from patho:  " + str(len(removed)))
    print("  Final harmful:       " + str(final_harmful))
    print("  Final uncertain:     " + str(final_uncertain))
    print("  Final benign:        " + str(final_benign))
    print("  Output: " + out_path)
    print("  Next: python code/step_16_clinvar_calibrate.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
