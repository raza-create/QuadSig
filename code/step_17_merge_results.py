#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 17: Merge All Signals Into results.csv")
    print("=" * 65)

    print("[Step 17] Loading all intermediate data...")

    struct = {}
    with open(os.path.join(DATA_INTER, "structural_scores.csv")) as f:
        for r in csv.DictReader(f): struct[r["mutation"]] = r

    pop = {}
    with open(os.path.join(DATA_INTER, "population_scores.csv")) as f:
        for r in csv.DictReader(f): pop[r["mutation"]] = r

    cons = {}
    with open(os.path.join(DATA_INTER, "conservation_scores.csv")) as f:
        for r in csv.DictReader(f): cons[int(r["position"])] = r

    phys = {}
    with open(os.path.join(DATA_INTER, "physicochemical_scores.csv")) as f:
        for r in csv.DictReader(f): phys[r["mutation"]] = r

    dom = {}
    with open(os.path.join(DATA_INTER, "domain_annotations.csv")) as f:
        for r in csv.DictReader(f): dom[r["mutation"]] = r

    zn = {}
    with open(os.path.join(DATA_INTER, "zinc_flags.csv")) as f:
        for r in csv.DictReader(f): zn[r["mutation"]] = r

    gof = {}
    with open(os.path.join(DATA_INTER, "gof_flags.csv")) as f:
        for r in csv.DictReader(f): gof[r["mutation"]] = r

    cal = {}
    with open(os.path.join(DATA_INTER, "calibrated_scores.csv")) as f:
        for r in csv.DictReader(f): cal[r["mutation"]] = r

    cv = {}
    cv_path = os.path.join(DATA_RAW, "clinvar_p53.csv")
    if os.path.exists(cv_path):
        with open(cv_path) as f:
            for r in csv.DictReader(f): cv[r["mutation"]] = r

    with open(os.path.join(DATA_INTER, "mutations.txt")) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 17] Total mutations: " + str(len(all_mutations)))

    results = []
    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        s = struct.get(mutation, {})
        p = pop.get(mutation, {})
        c = cons.get(pos, {})
        ph = phys.get(mutation, {})
        d = dom.get(mutation, {})
        z = zn.get(mutation, {})
        g = gof.get(mutation, {})
        ca = cal.get(mutation, {})
        cl = cv.get(mutation, {})

        results.append({
            "mutation": mutation,
            "position": pos,
            "wt_aa": mutation[0],
            "mut_aa": mutation[-1],
            "quadsig_score": ca.get("quadsig_score", ""),
            "category": ca.get("calibrated_category", ""),
            "population_score": p.get("population_score", ""),
            "allele_frequency": p.get("allele_frequency", ""),
            "conservation_score": c.get("conservation_score", ""),
            "rmsd": s.get("rmsd", ""),
            "plddt_wt": s.get("plddt_wt", ""),
            "plddt_mut": s.get("plddt_mut", ""),
            "plddt_change": s.get("plddt_change", ""),
            "physicochemical_score": ph.get("physicochemical_score", ""),
            "domain": d.get("domain_code", ""),
            "domain_name": d.get("domain_name", ""),
            "is_zinc_coord": z.get("is_zinc_coord", ""),
            "is_gof": g.get("is_gof", ""),
            "gof_mechanism": g.get("gof_mechanism", ""),
            "structural_divergence": g.get("structural_divergence_note", ""),
            "mechanism": ca.get("mechanism", ""),
            "clinvar_label": cl.get("clinvar_label", ""),
            "clinvar_significance": cl.get("clinvar_significance", ""),
            "clinvar_agreement": ca.get("clinvar_agreement", ""),
        })

    fields = list(results[0].keys())
    out_path = os.path.join(RESULTS_DIR, "results.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    lp = sum(1 for r in results if r["category"] == "likely_pathogenic")
    amb = sum(1 for r in results if r["category"] == "ambiguous")
    lb = sum(1 for r in results if r["category"] == "likely_benign")
    has_cv = sum(1 for r in results if r["clinvar_label"])
    zn_count = sum(1 for r in results if r["is_zinc_coord"] == "True")
    gof_count = sum(1 for r in results if r["is_gof"] == "True")

    all_scores = [float(r["quadsig_score"]) for r in results if r["quadsig_score"]]
    gene_avg = round(sum(all_scores) / len(all_scores), 4)

    print("[Step 17] Domain-level average QuadSig scores:")
    for code, (start, end, name) in DOMAINS.items():
        ds = [float(r["quadsig_score"]) for r in results if r["domain"] == code and r["quadsig_score"]]
        if ds:
            print("    " + code + ": " + str(round(sum(ds)/len(ds), 4)))

    print("")
    print("=" * 65)
    print("  STEP 17 COMPLETE — RESULTS.CSV GENERATED")
    print("=" * 65)
    print("  Rows:                " + str(len(results)))
    print("  Columns:             " + str(len(fields)))
    print("  File size:           " + str(os.path.getsize(out_path)) + " bytes")
    print("  Likely pathogenic:   " + str(lp))
    print("  Ambiguous:           " + str(amb))
    print("  Likely benign:       " + str(lb))
    print("  Gene-level score:    " + str(gene_avg))
    print("  ClinVar matches:     " + str(has_cv))
    print("  Zinc mutations:      " + str(zn_count))
    print("  GOF mutations:       " + str(gof_count))

    for hs in ["R175H","G245S","R248W","R248Q","R249S","R273H","R273C","P72R"]:
        for r in results:
            if r["mutation"] == hs:
                print("  " + hs + ": " + str(r["quadsig_score"]) + " " + r["category"] + " domain=" + r["domain"])
                break

    print("  Output: " + out_path)
    print("  Next: python code/step_18_validate.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
