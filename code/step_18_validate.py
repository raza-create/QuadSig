#!/usr/bin/env python3
import os, sys, csv, math

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def spearman(x, y):
    n = len(x)
    if n < 3:
        return 0
    def rank(lst):
        indexed = sorted(enumerate(lst), key=lambda t: t[1])
        ranks = [0] * n
        for rank_val, (orig_idx, val) in enumerate(indexed):
            ranks[orig_idx] = rank_val + 1
        return ranks
    rx = rank(x)
    ry = rank(y)
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
    return round(1 - 6 * d2 / (n * (n * n - 1)), 4)

def main():
    print("=" * 65)
    print("  QuadSig Step 18: Seven Independent Validation Benchmarks")
    print("=" * 65)

    # Load results
    results = {}
    with open(os.path.join(RESULTS_DIR, "results.csv")) as f:
        for r in csv.DictReader(f):
            results[r["mutation"]] = r
    print("[Step 18] Results loaded: " + str(len(results)))

    report = []

    # ── Benchmark 1: ClinVar Agreement ──
    print("\n[Benchmark 1] ClinVar Agreement")
    agree = 0
    total = 0
    for mut, r in results.items():
        if r["clinvar_agreement"] == "AGREE":
            agree += 1
            total += 1
        elif r["clinvar_agreement"] == "DISAGREE":
            total += 1
    cv_rate = round(agree / total * 100, 1) if total > 0 else 0
    target_met = cv_rate >= 85.0
    print("    Agreement: " + str(agree) + "/" + str(total) + " (" + str(cv_rate) + "%)")
    print("    Target: >= 85%  " + ("PASS" if target_met else "FAIL"))
    report.append({"benchmark": "ClinVar Agreement", "result": str(cv_rate) + "%", "target": ">= 85%", "status": "PASS" if target_met else "FAIL"})

    # ── Benchmark 2: Conservation Correlation ──
    print("\n[Benchmark 2] Conservation vs QuadSig Score Correlation")
    cons_scores = []
    quad_scores = []
    for mut, r in results.items():
        if r["conservation_score"] and r["quadsig_score"]:
            cons_scores.append(float(r["conservation_score"]))
            quad_scores.append(float(r["quadsig_score"]))
    cons_corr = spearman(cons_scores, quad_scores)
    target_met = cons_corr >= 0.50
    print("    Spearman correlation: " + str(cons_corr))
    print("    Target: >= 0.50  " + ("PASS" if target_met else "FAIL"))
    report.append({"benchmark": "Conservation Correlation", "result": str(cons_corr), "target": ">= 0.50", "status": "PASS" if target_met else "FAIL"})

    # ── Benchmark 3: Domain Gradient ──
    print("\n[Benchmark 3] Domain Score Gradient (DBD > tails)")
    domain_avgs = {}
    for code, (start, end, name) in DOMAINS.items():
        ds = [float(r["quadsig_score"]) for r in results.values() if r["domain"] == code and r["quadsig_score"]]
        if ds:
            domain_avgs[code] = round(sum(ds) / len(ds), 4)
    dbd = domain_avgs.get("DBD", 0)
    tad = domain_avgs.get("TAD", 0)
    ctd = domain_avgs.get("CTD", 0)
    prr = domain_avgs.get("PRR", 0)
    gradient = dbd > max(tad, ctd, prr)
    for code, avg in sorted(domain_avgs.items(), key=lambda x: -x[1]):
        print("    " + code + ": " + str(avg))
    print("    DBD highest: " + ("PASS" if gradient else "FAIL"))
    report.append({"benchmark": "Domain Gradient", "result": "DBD=" + str(dbd), "target": "DBD highest", "status": "PASS" if gradient else "FAIL"})

    # ── Benchmark 4: Zinc Residue Scores ──
    print("\n[Benchmark 4] Zinc Coordination Residue Scores")
    zinc_scores = [float(r["quadsig_score"]) for r in results.values() if r["is_zinc_coord"] == "True" and r["quadsig_score"]]
    zinc_avg = round(sum(zinc_scores) / len(zinc_scores), 4) if zinc_scores else 0
    zinc_min = round(min(zinc_scores), 4) if zinc_scores else 0
    zinc_above_90 = sum(1 for s in zinc_scores if s >= 0.90)
    zinc_above_63 = sum(1 for s in zinc_scores if s >= 0.63)
    target_met = zinc_avg >= 0.85
    print("    Zinc mutations: " + str(len(zinc_scores)))
    print("    Average score: " + str(zinc_avg))
    print("    Min score: " + str(zinc_min))
    print("    Above 0.90: " + str(zinc_above_90) + "/76")
    print("    Above 0.63 (pathogenic threshold): " + str(zinc_above_63) + "/76")
    print("    Target: avg >= 0.85  " + ("PASS" if target_met else "FAIL"))
    report.append({"benchmark": "Zinc Residue Scores", "result": "avg=" + str(zinc_avg), "target": "avg >= 0.85", "status": "PASS" if target_met else "FAIL"})

    # ── Benchmark 5: GOF Detection ──
    print("\n[Benchmark 5] Gain-of-Function Mutation Detection")
    gof_scores = []
    for gof_mut in GOF_MUTATIONS:
        r = results.get(gof_mut)
        if r:
            s = float(r["quadsig_score"])
            gof_scores.append(s)
            print("    " + gof_mut + ": " + str(s) + " " + r["category"])
    all_patho = all(s >= 0.63 for s in gof_scores)
    print("    All GOF pathogenic: " + ("PASS" if all_patho else "FAIL"))
    report.append({"benchmark": "GOF Detection", "result": "All >= 0.63: " + str(all_patho), "target": "All pathogenic", "status": "PASS" if all_patho else "FAIL"})

    # ── Benchmark 6: BLOSUM62 Anti-correlation ──
    print("\n[Benchmark 6] BLOSUM62 Anti-correlation")
    BLOSUM62 = {('A','A'):4,('R','R'):5,('N','N'):6,('D','D'):6,('C','C'):9,('Q','Q'):5,('E','E'):5,('G','G'):6,('H','H'):8,('I','I'):4,('L','L'):4,('K','K'):5,('M','M'):5,('F','F'):6,('P','P'):7,('S','S'):4,('T','T'):5,('W','W'):11,('Y','Y'):7,('V','V'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,('A','Q'):-1,('A','E'):-1,('A','G'):0,('A','H'):-2,('A','I'):-1,('A','L'):-1,('A','K'):-1,('A','M'):-1,('A','F'):-2,('A','P'):-1,('A','S'):1,('A','T'):0,('A','W'):-3,('A','Y'):-2,('A','V'):0,('R','N'):0,('R','D'):-2,('R','C'):-3,('R','Q'):1,('R','E'):0,('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,('R','K'):2,('R','M'):-1,('R','F'):-3,('R','P'):-2,('R','S'):-1,('R','T'):-1,('R','W'):-3,('R','Y'):-2,('R','V'):-3,('N','D'):1,('N','C'):-3,('N','Q'):0,('N','E'):0,('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,('N','K'):0,('N','M'):-2,('N','F'):-3,('N','P'):-2,('N','S'):1,('N','T'):0,('N','W'):-4,('N','Y'):-2,('N','V'):-3,('D','C'):-3,('D','Q'):0,('D','E'):2,('D','G'):-1,('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,('D','M'):-3,('D','F'):-3,('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,('D','Y'):-3,('D','V'):-3,('C','Q'):-3,('C','E'):-4,('C','G'):-3,('C','H'):-3,('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,('Q','L'):-2,('Q','K'):1,('Q','M'):0,('Q','F'):-3,('Q','P'):-1,('Q','S'):0,('Q','T'):-1,('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,('E','G'):-2,('E','H'):0,('E','I'):-3,('E','L'):-3,('E','K'):1,('E','M'):-2,('E','F'):-3,('E','P'):-1,('E','S'):0,('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,('G','M'):-3,('G','F'):-3,('G','P'):-2,('G','S'):0,('G','T'):-2,('G','W'):-2,('G','Y'):-3,('G','V'):-3,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,('H','F'):-1,('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,('H','Y'):2,('H','V'):-3,('I','L'):2,('I','K'):-3,('I','M'):1,('I','F'):0,('I','P'):-3,('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,('L','K'):-2,('L','M'):2,('L','F'):0,('L','P'):-3,('L','S'):-2,('L','T'):-1,('L','W'):-2,('L','Y'):-1,('L','V'):1,('K','M'):-1,('K','F'):-3,('K','P'):-1,('K','S'):0,('K','T'):-1,('K','W'):-3,('K','Y'):-2,('K','V'):-2,('M','F'):0,('M','P'):-2,('M','S'):-1,('M','T'):-1,('M','W'):-1,('M','Y'):-1,('M','V'):1,('F','P'):-4,('F','S'):-2,('F','T'):-2,('F','W'):1,('F','Y'):3,('F','V'):-1,('P','S'):-1,('P','T'):-1,('P','W'):-4,('P','Y'):-3,('P','V'):-2,('S','T'):1,('S','W'):-3,('S','Y'):-2,('S','V'):-2,('T','W'):-3,('T','Y'):-2,('T','V'):0,('W','Y'):2,('W','V'):-3,('Y','V'):-1}
    def blosum_score(a, b):
        if (a, b) in BLOSUM62: return BLOSUM62[(a, b)]
        if (b, a) in BLOSUM62: return BLOSUM62[(b, a)]
        return -1
    bl_scores = []
    qs_scores = []
    for mut, r in results.items():
        if r["quadsig_score"]:
            bl = blosum_score(r["wt_aa"], r["mut_aa"])
            bl_scores.append(bl)
            qs_scores.append(float(r["quadsig_score"]))
    bl_corr = spearman(bl_scores, qs_scores)
    target_met = bl_corr <= -0.30
    print("    Spearman: " + str(bl_corr))
    print("    Target: <= -0.30  " + ("PASS" if target_met else "FAIL"))
    report.append({"benchmark": "BLOSUM62 Anti-correlation", "result": str(bl_corr), "target": "<= -0.30", "status": "PASS" if target_met else "FAIL"})

    # ── Benchmark 7: Pathogenic Enrichment in DBD ──
    print("\n[Benchmark 7] Pathogenic Enrichment in DBD")
    dbd_patho = sum(1 for r in results.values() if r["domain"] == "DBD" and r["category"] == "likely_pathogenic")
    dbd_total = sum(1 for r in results.values() if r["domain"] == "DBD")
    dbd_rate = round(dbd_patho / dbd_total * 100, 1) if dbd_total > 0 else 0
    non_dbd_patho = sum(1 for r in results.values() if r["domain"] != "DBD" and r["domain"] != "" and r["category"] == "likely_pathogenic")
    non_dbd_total = sum(1 for r in results.values() if r["domain"] != "DBD" and r["domain"] != "")
    non_dbd_rate = round(non_dbd_patho / non_dbd_total * 100, 1) if non_dbd_total > 0 else 0
    enrichment = dbd_rate > non_dbd_rate
    print("    DBD pathogenic rate: " + str(dbd_rate) + "%")
    print("    Non-DBD pathogenic rate: " + str(non_dbd_rate) + "%")
    print("    DBD enriched: " + ("PASS" if enrichment else "FAIL"))
    report.append({"benchmark": "DBD Pathogenic Enrichment", "result": "DBD=" + str(dbd_rate) + "% vs " + str(non_dbd_rate) + "%", "target": "DBD > others", "status": "PASS" if enrichment else "FAIL"})

    # ── Save validation report ──
    report_path = os.path.join(RESULTS_DIR, "validation_report.csv")
    with open(report_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["benchmark","result","target","status"])
        w.writeheader()
        w.writerows(report)

    passed = sum(1 for r in report if r["status"] == "PASS")
    total_b = len(report)

    print("")
    print("=" * 65)
    print("  STEP 18 COMPLETE — VALIDATION REPORT")
    print("=" * 65)
    for r in report:
        print("    " + r["status"] + "  " + r["benchmark"] + ": " + r["result"] + " (target: " + r["target"] + ")")
    print("")
    print("  PASSED: " + str(passed) + " / " + str(total_b))
    print("  Output: " + report_path)
    print("  Next: python code/step_19_figures.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
