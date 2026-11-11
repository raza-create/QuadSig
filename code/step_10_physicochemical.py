#!/usr/bin/env python3
import os, sys, csv, math

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 10: Physicochemical Disruption Scores")
    print("=" * 65)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 10] Total mutations: " + str(len(all_mutations)))

    # Compute max possible differences for normalization
    charges = list(AA_CHARGE.values())
    mws = list(AA_MW.values())
    hydros = list(AA_HYDRO.values())
    max_charge = max(charges) - min(charges)
    max_mw = max(mws) - min(mws)
    max_hydro = max(hydros) - min(hydros)

    print("[Step 10] Normalization ranges:")
    print("    Charge: " + str(max_charge))
    print("    Size (MW): " + str(max_mw))
    print("    Hydrophobicity: " + str(max_hydro))

    results = []
    for mutation in all_mutations:
        wt_aa = mutation[0]
        mut_aa = mutation[-1]
        pos = int(mutation[1:-1])

        charge_diff = abs(AA_CHARGE[mut_aa] - AA_CHARGE[wt_aa]) / max_charge if max_charge > 0 else 0
        size_diff = abs(AA_MW[mut_aa] - AA_MW[wt_aa]) / max_mw if max_mw > 0 else 0
        hydro_diff = abs(AA_HYDRO[mut_aa] - AA_HYDRO[wt_aa]) / max_hydro if max_hydro > 0 else 0

        # Combined score: equal weight, Euclidean distance normalized to 0-1
        combined = math.sqrt((charge_diff**2 + size_diff**2 + hydro_diff**2) / 3)
        combined = round(min(combined, 1.0), 4)

        results.append({
            "mutation": mutation,
            "position": pos,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "charge_diff": round(charge_diff, 4),
            "size_diff": round(size_diff, 4),
            "hydrophobicity_diff": round(hydro_diff, 4),
            "physicochemical_score": combined,
        })

    out_path = os.path.join(DATA_INTER, "physicochemical_scores.csv")
    fields = ["mutation","position","wt_aa","mut_aa","charge_diff","size_diff","hydrophobicity_diff","physicochemical_score"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    scores = [r["physicochemical_score"] for r in results]
    avg = round(sum(scores) / len(scores), 4)

    print("")
    print("[Step 10] Example scores (from proposal):")
    for hs in ["R175H","R248W","R248Q","R273H","R273C","G245S"]:
        for r in results:
            if r["mutation"] == hs:
                print("    " + hs + ": charge=" + str(r["charge_diff"]) + " size=" + str(r["size_diff"]) + " hydro=" + str(r["hydrophobicity_diff"]) + " combined=" + str(r["physicochemical_score"]))
                break

    # Proposal says R->W ~0.91, R->K ~0.21
    for test in ["R175W","R175K"]:
        for r in results:
            if r["mutation"] == test:
                print("    " + test + ": " + str(r["physicochemical_score"]) + " (proposal: R->W ~0.91, R->K ~0.21)")
                break

    print("")
    print("=" * 65)
    print("  STEP 10 COMPLETE")
    print("=" * 65)
    print("  Mutations scored:  " + str(len(results)))
    print("  Average score:     " + str(avg))
    print("  Min score:         " + str(round(min(scores), 4)))
    print("  Max score:         " + str(round(max(scores), 4)))
    print("  Output: " + out_path + " (" + str(os.path.getsize(out_path)) + " bytes)")
    print("")
    print("  Next: python code/step_11_domain_annotation.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
