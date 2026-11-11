#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 13: Gain-of-Function (GOF) Mutation Detection")
    print("=" * 65)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 13] Total mutations: " + str(len(all_mutations)))
    print("[Step 13] Known GOF hotspots: " + str(len(GOF_MUTATIONS)))
    for gof, mech in GOF_MUTATIONS.items():
        print("    " + gof + ": " + mech)

    # Load structural and physicochemical scores for divergence analysis
    struct_data = {}
    struct_path = os.path.join(DATA_INTER, "structural_scores.csv")
    if os.path.exists(struct_path):
        with open(struct_path) as f:
            for row in csv.DictReader(f):
                struct_data[row["mutation"]] = row

    physchem_data = {}
    physchem_path = os.path.join(DATA_INTER, "physicochemical_scores.csv")
    if os.path.exists(physchem_path):
        with open(physchem_path) as f:
            for row in csv.DictReader(f):
                physchem_data[row["mutation"]] = row

    results = []
    gof_count = 0
    for mutation in all_mutations:
        is_gof = mutation in GOF_MUTATIONS
        mechanism = GOF_MUTATIONS.get(mutation, "")

        # Structural divergence note for GOF mutations
        divergence = ""
        if is_gof:
            gof_count += 1
            sd = struct_data.get(mutation, {})
            rmsd = float(sd.get("rmsd", 0))
            pc = physchem_data.get(mutation, {})
            pc_score = float(pc.get("physicochemical_score", 0))

            if rmsd < 1.0:
                divergence = "Low RMSD despite pathogenicity — GOF acts via surface chemistry change, not structural disruption"
            elif rmsd < 2.0:
                divergence = "Moderate RMSD — partial structural change plus GOF surface interaction"
            else:
                divergence = "High RMSD — structural disruption plus GOF gain of novel interactions"

        results.append({
            "mutation": mutation,
            "position": int(mutation[1:-1]),
            "wt_aa": mutation[0],
            "mut_aa": mutation[-1],
            "is_gof": is_gof,
            "gof_mechanism": mechanism,
            "structural_divergence_note": divergence,
        })

    out_path = os.path.join(DATA_INTER, "gof_flags.csv")
    fields = ["mutation","position","wt_aa","mut_aa","is_gof","gof_mechanism","structural_divergence_note"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    print("")
    print("[Step 13] GOF mutation structural analysis:")
    for gof_mut in GOF_MUTATIONS.keys():
        sd = struct_data.get(gof_mut, {})
        pc = physchem_data.get(gof_mut, {})
        rmsd = sd.get("rmsd", "?")
        plddt_ch = sd.get("plddt_change", "?")
        pc_score = pc.get("physicochemical_score", "?")
        for r in results:
            if r["mutation"] == gof_mut:
                print("    " + gof_mut + ": RMSD=" + str(rmsd) + " pLDDT_change=" + str(plddt_ch) + " physchem=" + str(pc_score))
                print("      " + r["structural_divergence_note"])
                break

    print("")
    print("=" * 65)
    print("  STEP 13 COMPLETE")
    print("=" * 65)
    print("  GOF mutations flagged: " + str(gof_count) + " / 7")
    print("  Non-GOF mutations:     " + str(len(results) - gof_count))
    print("  Output: " + out_path)
    print("  Next: python code/step_14_fusion_scoring.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
