#!/usr/bin/env python3
import os, sys, csv, requests, time

sys.path.insert(0, os.path.dirname(__file__))
from config import *

AA_3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}

def parse_pchange(title):
    if "(p." not in title:
        return ""
    pc = title.split("(p.")[1].split(")")[0]
    if len(pc) < 6:
        return ""
    wt3 = pc[:3]
    mut3 = pc[-3:]
    pos_str = pc[3:-3]
    wt = AA_3TO1.get(wt3, "")
    mut = AA_3TO1.get(mut3, "")
    if wt and mut and pos_str.isdigit():
        return wt + pos_str + mut
    return ""

def main():
    print("=" * 65)
    print("  QuadSig Step 16: Calibrate Scores Using ClinVar")
    print("=" * 65)

    scores = {}
    with open(os.path.join(DATA_INTER, "quadsig_scores.csv")) as f:
        for r in csv.DictReader(f):
            scores[r["mutation"]] = float(r["quadsig_score"])
    print("[Step 16] QuadSig scores loaded: " + str(len(scores)))

    # Search ClinVar for TP53 variants
    print("[Step 16] Searching ClinVar for TP53 variants...")
    url = CLINVAR_API + "/esearch.fcgi?db=clinvar&term=TP53[gene]+AND+single_nucleotide_variant[Type+of+variation]&retmax=3000&retmode=json"
    resp = requests.get(url, timeout=60)
    data = resp.json()
    ids = data.get("esearchresult", {}).get("idlist", [])
    print("[Step 16] ClinVar IDs found: " + str(len(ids)))

    # Fetch details in batches
    print("[Step 16] Fetching variant details...")
    clinvar_variants = []
    for i in range(0, len(ids), 100):
        batch = ids[i:i+100]
        furl = CLINVAR_API + "/esummary.fcgi?db=clinvar&id=" + ",".join(batch) + "&retmode=json"
        try:
            resp = requests.get(furl, timeout=60)
            result = resp.json().get("result", {})
            for vid in batch:
                entry = result.get(vid, {})
                title = entry.get("title", "")
                gene_list = entry.get("genes", [])
                gene = gene_list[0].get("symbol", "") if gene_list else ""
                if gene != "TP53":
                    continue

                sig = entry.get("germline_classification", {}).get("description", "")
                if not sig:
                    continue

                mutation = parse_pchange(title)
                if not mutation:
                    continue

                sig_lower = sig.lower()
                if "pathogenic" in sig_lower and "benign" not in sig_lower and "uncertain" not in sig_lower and "conflicting" not in sig_lower:
                    label = "pathogenic"
                elif "likely pathogenic" in sig_lower:
                    label = "pathogenic"
                elif "benign" in sig_lower and "pathogenic" not in sig_lower:
                    label = "benign"
                elif "likely benign" in sig_lower:
                    label = "benign"
                else:
                    label = "uncertain"

                clinvar_variants.append({
                    "mutation": mutation,
                    "clinvar_significance": sig,
                    "clinvar_label": label,
                })
        except Exception as e:
            print("    [ERROR] batch " + str(i) + ": " + str(e))

        if (i + 100) % 500 == 0:
            print("    Fetched " + str(min(i + 100, len(ids))) + "/" + str(len(ids)))
        time.sleep(0.5)

    # Deduplicate by mutation
    clinvar_dict = {}
    for cv in clinvar_variants:
        mut = cv["mutation"]
        if mut not in clinvar_dict:
            clinvar_dict[mut] = cv

    print("[Step 16] Unique ClinVar variants: " + str(len(clinvar_dict)))

    # Save ClinVar data
    clinvar_path = os.path.join(DATA_RAW, "clinvar_p53.csv")
    with open(clinvar_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mutation","clinvar_significance","clinvar_label"])
        w.writeheader()
        for cv in clinvar_dict.values():
            w.writerow(cv)
    print("[Step 16] Saved: " + clinvar_path)

    # Match with our scores
    matched_patho = []
    matched_benign = []
    for mut, cv in clinvar_dict.items():
        if mut in scores:
            if cv["clinvar_label"] == "pathogenic":
                matched_patho.append(scores[mut])
            elif cv["clinvar_label"] == "benign":
                matched_benign.append(scores[mut])

    print("[Step 16] Matched pathogenic: " + str(len(matched_patho)))
    print("[Step 16] Matched benign: " + str(len(matched_benign)))

    # Calibrate thresholds
    patho_threshold = PATHOGENIC_THRESHOLD
    benign_threshold = BENIGN_THRESHOLD

    if matched_patho:
        matched_patho.sort()
        p10 = matched_patho[max(0, int(len(matched_patho) * 0.10))]
        print("[Step 16] Pathogenic 10th percentile: " + str(round(p10, 4)))
        patho_threshold = round(max(p10 - 0.02, 0.50), 2)

    if matched_benign:
        matched_benign.sort()
        b90 = matched_benign[min(len(matched_benign)-1, int(len(matched_benign) * 0.90))]
        print("[Step 16] Benign 90th percentile: " + str(round(b90, 4)))
        benign_threshold = round(min(b90 + 0.02, 0.50), 2)

    if patho_threshold <= benign_threshold:
        patho_threshold = 0.66
        benign_threshold = 0.34

    print("[Step 16] Calibrated thresholds:")
    print("    Likely pathogenic: >= " + str(patho_threshold))
    print("    Ambiguous: " + str(benign_threshold) + " - " + str(patho_threshold))
    print("    Likely benign: < " + str(benign_threshold))

    # Apply to all mutations
    results = []
    with open(os.path.join(DATA_INTER, "quadsig_scores.csv")) as f:
        for r in csv.DictReader(f):
            mut = r["mutation"]
            score = float(r["quadsig_score"])
            if score >= patho_threshold:
                cat = "likely_pathogenic"
            elif score < benign_threshold:
                cat = "likely_benign"
            else:
                cat = "ambiguous"

            cv = clinvar_dict.get(mut, {})
            cv_label = cv.get("clinvar_label", "")
            agreement = ""
            if cv_label == "pathogenic" and cat == "likely_pathogenic":
                agreement = "AGREE"
            elif cv_label == "benign" and cat == "likely_benign":
                agreement = "AGREE"
            elif cv_label in ["pathogenic","benign"]:
                agreement = "DISAGREE"

            results.append({
                "mutation": mut,
                "quadsig_score": score,
                "calibrated_category": cat,
                "clinvar_label": cv_label,
                "clinvar_agreement": agreement,
                "mechanism": r.get("mechanism", ""),
            })

    out_path = os.path.join(DATA_INTER, "calibrated_scores.csv")
    fields = ["mutation","quadsig_score","calibrated_category","clinvar_label","clinvar_agreement","mechanism"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    lp = sum(1 for r in results if r["calibrated_category"] == "likely_pathogenic")
    amb = sum(1 for r in results if r["calibrated_category"] == "ambiguous")
    lb = sum(1 for r in results if r["calibrated_category"] == "likely_benign")
    agree = sum(1 for r in results if r["clinvar_agreement"] == "AGREE")
    disagree = sum(1 for r in results if r["clinvar_agreement"] == "DISAGREE")
    total_j = agree + disagree
    rate = round(agree / total_j * 100, 1) if total_j > 0 else 0

    print("")
    print("=" * 65)
    print("  STEP 16 COMPLETE")
    print("=" * 65)
    print("  Likely pathogenic: " + str(lp))
    print("  Ambiguous:         " + str(amb))
    print("  Likely benign:     " + str(lb))
    print("  ClinVar agreement: " + str(agree) + "/" + str(total_j) + " (" + str(rate) + "%) target > 85%")

    for hs in ["R175H","G245S","R248W","R273H","P72R"]:
        for r in results:
            if r["mutation"] == hs:
                print("  " + hs + ": " + str(r["quadsig_score"]) + " -> " + r["calibrated_category"] + " ClinVar=" + r["clinvar_label"] + " " + r["clinvar_agreement"])
                break

    print("  Output: " + out_path)
    print("  Next: python code/step_17_merge_results.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
