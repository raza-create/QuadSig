#!/usr/bin/env python3
import os, sys, csv, math

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 14: Fuse All Signals - Enhanced Formula Scoring")
    print("=" * 65)

    # Load all signals
    struct = {}
    with open(os.path.join(DATA_INTER, "structural_scores.csv")) as f:
        for r in csv.DictReader(f):
            struct[r["mutation"]] = r
    pop = {}
    with open(os.path.join(DATA_INTER, "population_scores.csv")) as f:
        for r in csv.DictReader(f):
            pop[r["mutation"]] = r
    cons = {}
    with open(os.path.join(DATA_INTER, "conservation_scores.csv")) as f:
        for r in csv.DictReader(f):
            cons[int(r["position"])] = float(r["conservation_score"])
    phys = {}
    with open(os.path.join(DATA_INTER, "physicochemical_scores.csv")) as f:
        for r in csv.DictReader(f):
            phys[r["mutation"]] = r
    domain = {}
    with open(os.path.join(DATA_INTER, "domain_annotations.csv")) as f:
        for r in csv.DictReader(f):
            domain[r["mutation"]] = r
    zinc = {}
    with open(os.path.join(DATA_INTER, "zinc_flags.csv")) as f:
        for r in csv.DictReader(f):
            zinc[r["mutation"]] = r
    gof = {}
    with open(os.path.join(DATA_INTER, "gof_flags.csv")) as f:
        for r in csv.DictReader(f):
            gof[r["mutation"]] = r

    with open(os.path.join(DATA_INTER, "mutations.txt")) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 14] Total mutations: " + str(len(all_mutations)))

    results = []
    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        s = struct.get(mutation, {})
        p = pop.get(mutation, {})
        ph = phys.get(mutation, {})
        d = domain.get(mutation, {})
        z = zinc.get(mutation, {})
        g = gof.get(mutation, {})

        pop_score = float(p.get("population_score", 1.0))
        cons_score = cons.get(pos, 0.5)
        rmsd = float(s.get("rmsd", 0))
        plddt_change = float(s.get("plddt_change", 0))
        phys_score = float(ph.get("physicochemical_score", 0))
        dom = d.get("domain_code", "LINKER")
        is_zinc = z.get("is_zinc_coord", "False") == "True"
        is_gof = g.get("is_gof", "False") == "True"

        # Normalize RMSD: sigmoid-like mapping
        # Small RMSD (<0.5) = low signal, large RMSD (>3) = strong signal
        rmsd_norm = 1.0 / (1.0 + math.exp(-1.5 * (rmsd - 2.0)))

        # Normalize pLDDT change: negative change = bad
        plddt_penalty = min(abs(plddt_change) / 30.0, 1.0) if plddt_change < 0 else 0

        # Domain weight: DBD mutations more likely pathogenic
        domain_weight = {"DBD": 0.15, "TET": 0.08, "TAD": 0.02, "PRR": 0.0, "CTD": -0.05, "LINKER": -0.03}.get(dom, 0)

        # Conservation interaction: high conservation + high disruption = very pathogenic
        cons_phys_interaction = cons_score * phys_score * 0.5

        # Base score: weighted combination
        base = (
            0.15 * pop_score +
            0.25 * cons_score +
            0.20 * rmsd_norm +
            0.10 * plddt_penalty +
            0.15 * phys_score +
            0.15 * cons_phys_interaction
        )

        # Add domain context
        score = base + domain_weight

        # Population override: if variant is common (AF > 0.001), strongly push toward benign
        af = float(p.get("allele_frequency", 0))
        if af > 0.01:
            score = score * 0.2
        elif af > 0.001:
            score = score * 0.4
        elif af > 0.0001:
            score = score * 0.6

        # Zinc coordination: always pathogenic
        if is_zinc:
            score = max(score, 0.93)

        # GOF: always pathogenic
        if is_gof:
            score = max(score, 0.91)

        # Low conservation in disordered region: push toward benign
        if cons_score < 0.3 and dom in ["CTD", "TAD", "PRR", "LINKER"] and rmsd < 1.5:
            score = score * 0.5

        # High conservation + DBD + any structural change: push toward pathogenic
        if cons_score > 0.8 and dom == "DBD" and rmsd > 0.5:
            score = max(score, score + 0.1)

        score = round(max(0.0, min(1.0, score)), 4)

        if score >= 0.66:
            verdict = "harmful"
        elif score >= 0.34:
            verdict = "uncertain"
        else:
            verdict = "benign"

        # Generate mechanism text
        mech_parts = []
        if is_zinc:
            mech_parts.append("Zinc coordination disrupted at " + mutation[0] + str(pos))
        if is_gof:
            mech_parts.append("Known GOF mutation: " + g.get("gof_mechanism", ""))
        if dom == "DBD":
            mech_parts.append("DNA-binding domain")
        if cons_score > 0.8:
            mech_parts.append("Highly conserved position (cons=" + str(round(cons_score, 2)) + ")")
        if rmsd > 2.0:
            mech_parts.append("Major structural change (RMSD=" + str(round(rmsd, 2)) + ")")
        if plddt_change < -10:
            mech_parts.append("Confidence drop (pLDDT change=" + str(round(plddt_change, 1)) + ")")
        if phys_score > 0.5:
            mech_parts.append("Chemically disruptive substitution")
        if af > 0.01:
            mech_parts.append("Common polymorphism (AF=" + str(round(af, 4)) + ")")
        if cons_score < 0.3:
            mech_parts.append("Variable position across species")
        if dom in ["CTD", "TAD"] and cons_score < 0.4:
            mech_parts.append("Disordered tail region")

        mechanism = "; ".join(mech_parts) if mech_parts else "No strong signals"

        results.append({
            "mutation": mutation,
            "quadsig_score": score,
            "verdict": verdict,
            "mechanism": mechanism,
        })

    out_path = os.path.join(DATA_INTER, "quadsig_scores.csv")
    fields = ["mutation", "quadsig_score", "verdict", "mechanism"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    scores_list = [r["quadsig_score"] for r in results]
    harmful = sum(1 for r in results if r["verdict"] == "harmful")
    uncertain = sum(1 for r in results if r["verdict"] == "uncertain")
    benign = sum(1 for r in results if r["verdict"] == "benign")

    print("")
    print("=" * 65)
    print("  STEP 14 COMPLETE (Enhanced Formula)")
    print("=" * 65)
    print("  Scored:     " + str(len(results)))
    print("  Harmful:    " + str(harmful))
    print("  Uncertain:  " + str(uncertain))
    print("  Benign:     " + str(benign))
    print("  Avg score:  " + str(round(sum(scores_list)/len(scores_list), 4)))

    for hs in ["R175H","G245S","R248W","R273H","R273C","P72R"]:
        for r in results:
            if r["mutation"] == hs:
                print("  " + hs + ": " + str(r["quadsig_score"]) + " (" + r["verdict"] + ") -- " + r["mechanism"][:80])
                break

    print("  Output: " + out_path)
    print("  Next: python code/step_15_self_distillation.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
