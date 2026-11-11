#!/usr/bin/env python3
import os, sys, csv, time, requests, json

sys.path.insert(0, os.path.dirname(__file__))
from config import *

CODON_TO_AA = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

COMP = {'A':'T','T':'A','C':'G','G':'C'}

def main():
    print("=" * 65)
    print("  QuadSig Step 7: Query gnomAD for Allele Frequencies")
    print("=" * 65)

    coords_path = os.path.join(DATA_INTER, "genomic_coords.csv")
    with open(coords_path) as f:
        coords = list(csv.DictReader(f))
    print("[Step 7] Total mutations: " + str(len(coords)))

    # Get CDS sequence
    print("[Step 7] Fetching CDS sequence...")
    resp = requests.get(ENSEMBL_API + "/sequence/id/ENST00000269305?type=cds", headers={"Content-Type":"application/json"}, timeout=60)
    resp.raise_for_status()
    cds = resp.json()["seq"]
    print("[Step 7] CDS: " + str(len(cds)) + " bp")

    # Get genomic reference sequence
    print("[Step 7] Fetching genomic reference...")
    all_gpos = set()
    for row in coords:
        if row["genomic_start"] != "NA":
            gs = int(row["genomic_start"])
            ge = int(row["genomic_end"])
            for p in range(gs, ge + 1):
                all_gpos.add(p)
    gmin = min(all_gpos)
    gmax = max(all_gpos)

    resp = requests.get(ENSEMBL_API + "/sequence/region/human/17:" + str(gmin) + ".." + str(gmax) + ":1", headers={"Content-Type":"application/json"}, timeout=60)
    resp.raise_for_status()
    genomic_seq = resp.json()["seq"]
    print("[Step 7] Genomic range: " + str(gmin) + "-" + str(gmax) + " (" + str(len(genomic_seq)) + " bp)")

    def get_ref(gpos):
        idx = gpos - gmin
        if 0 <= idx < len(genomic_seq):
            return genomic_seq[idx]
        return None

    print("[Step 7] Querying gnomAD (correct negative strand codon mapping)...")

    results = []
    query_count = 0
    found_count = 0

    for i, row in enumerate(coords):
        mutation = row["mutation"]
        pos = int(row["position"])
        wt_aa = row["wt_aa"]
        mut_aa = row["mut_aa"]
        gs_str = row["genomic_start"]

        if gs_str == "NA":
            results.append({"mutation": mutation, "position": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "allele_frequency": 0, "population_score": 1.0})
            continue

        gs = int(gs_str)
        ge = int(row["genomic_end"])

        # NEGATIVE STRAND: codon pos 0 = highest genomic, pos 2 = lowest
        codon_gpos = [ge, ge - 1, gs]  # codon order: 1st, 2nd, 3rd base

        ref_bases = [get_ref(p) for p in codon_gpos]
        if None in ref_bases:
            results.append({"mutation": mutation, "position": pos, "wt_aa": wt_aa, "mut_aa": mut_aa, "allele_frequency": 0, "population_score": 1.0})
            continue

        # CDS codon = complement of each genomic base (read in codon order)
        cds_codon = "".join(COMP[b] for b in ref_bases)
        encoded_aa = CODON_TO_AA.get(cds_codon, "?")

        best_af = 0.0

        # Try single nucleotide changes at each codon position
        for ci in range(3):
            gp = codon_gpos[ci]
            ref_g = ref_bases[ci]

            for alt_g in ['A', 'C', 'G', 'T']:
                if alt_g == ref_g:
                    continue

                # Make mutant CDS codon
                mut_bases = list(ref_bases)
                mut_bases[ci] = alt_g
                mut_cds_codon = "".join(COMP[b] for b in mut_bases)
                result_aa = CODON_TO_AA.get(mut_cds_codon, "?")

                if result_aa == mut_aa:
                    vid = "17-" + str(gp) + "-" + ref_g + "-" + alt_g
                    try:
                        q = '{variant(variantId: "' + vid + '", dataset: gnomad_r4) { exome { af } genome { af }}}'
                        resp = requests.post(GNOMAD_API, json={"query": q}, timeout=15)
                        if resp.status_code == 200:
                            data = resp.json()
                            v = data.get("data", {}).get("variant")
                            if v:
                                eaf = (v.get("exome") or {}).get("af", 0) or 0
                                gaf = (v.get("genome") or {}).get("af", 0) or 0
                                af = max(eaf, gaf)
                                if af > best_af:
                                    best_af = af
                        query_count += 1
                        time.sleep(0.05)
                    except:
                        query_count += 1
                        continue

        if best_af > 0:
            found_count += 1

        if best_af == 0:
            pop_score = 1.0
        elif best_af < 0.00001:
            pop_score = 0.95
        elif best_af < 0.0001:
            pop_score = 0.8
        elif best_af < 0.001:
            pop_score = 0.5
        elif best_af < 0.01:
            pop_score = 0.2
        else:
            pop_score = 0.0

        results.append({
            "mutation": mutation,
            "position": pos,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "allele_frequency": best_af,
            "population_score": round(pop_score, 4),
        })

        if (i + 1) % 500 == 0:
            print("    " + str(i + 1) + "/" + str(len(coords)) + " | queries: " + str(query_count) + " | found: " + str(found_count))

    out_path = os.path.join(DATA_INTER, "population_scores.csv")
    fields = ["mutation","position","wt_aa","mut_aa","allele_frequency","population_score"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    zero_af = sum(1 for r in results if r["allele_frequency"] == 0)
    nonzero = sum(1 for r in results if r["allele_frequency"] > 0)

    print("")
    print("=" * 65)
    print("  STEP 7 COMPLETE")
    print("=" * 65)
    print("  Total scored:      " + str(len(results)))
    print("  gnomAD queries:    " + str(query_count))
    print("  AF = 0 (unseen):   " + str(zero_af))
    print("  AF > 0 (observed): " + str(nonzero))

    for hs in ["R175H","G245S","R248W","R273H","P72R","P72A"]:
        for r in results:
            if r["mutation"] == hs:
                print("  " + hs + ": AF=" + str(r["allele_frequency"]) + " score=" + str(r["population_score"]))
                break

    print("")
    print("  Output: " + out_path)
    print("  Next: python code/step_08_download_msa.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
