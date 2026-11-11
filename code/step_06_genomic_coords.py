#!/usr/bin/env python3
import os, sys, csv, time, requests, json

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 6: Protein to Genomic Coordinates")
    print("=" * 65)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 6] Total mutations: " + str(len(all_mutations)))

    print("[Step 6] Fetching TP53 gene info from Ensembl...")
    url = ENSEMBL_API + "/lookup/symbol/homo_sapiens/TP53?expand=1"
    resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
    resp.raise_for_status()
    gene = resp.json()
    chrom = gene.get("seq_region_name", "17")
    strand = gene.get("strand", -1)
    canonical = gene.get("canonical_transcript", "").split(".")[0]
    print("[Step 6] Transcript: " + canonical + " chr" + chrom + " strand=" + str(strand))

    # Map each protein position via CDS coordinates
    # Protein pos N -> CDS positions (N-1)*3+1 to N*3
    print("[Step 6] Mapping 393 protein positions via CDS->genomic...")

    pos_to_genomic = {}
    for pos in range(1, SEQUENCE_LEN + 1):
        cds_start = (pos - 1) * 3 + 1
        cds_end = pos * 3
        url = ENSEMBL_API + "/map/cds/" + canonical + "/" + str(cds_start) + ".." + str(cds_end)
        try:
            resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            mappings = data.get("mappings", [])
            if mappings:
                gstarts = [m["start"] for m in mappings]
                gends = [m["end"] for m in mappings]
                pos_to_genomic[pos] = {
                    "start": min(gstarts),
                    "end": max(gends),
                    "chrom": mappings[0].get("seq_region_name", chrom),
                }
            else:
                pos_to_genomic[pos] = None
        except Exception as e:
            pos_to_genomic[pos] = None

        if pos % 50 == 0:
            print("    " + str(pos) + "/393 mapped")
        time.sleep(0.1)

    mapped_count = sum(1 for v in pos_to_genomic.values() if v is not None)
    print("[Step 6] Mapped " + str(mapped_count) + "/393 positions")

    results = []
    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        ginfo = pos_to_genomic.get(pos)
        if ginfo:
            results.append({
                "mutation": mutation,
                "position": pos,
                "wt_aa": mutation[0],
                "mut_aa": mutation[-1],
                "chromosome": ginfo["chrom"],
                "genomic_start": ginfo["start"],
                "genomic_end": ginfo["end"],
                "strand": strand,
            })
        else:
            results.append({
                "mutation": mutation,
                "position": pos,
                "wt_aa": mutation[0],
                "mut_aa": mutation[-1],
                "chromosome": chrom,
                "genomic_start": "NA",
                "genomic_end": "NA",
                "strand": strand,
            })

    out_path = os.path.join(DATA_INTER, "genomic_coords.csv")
    fields = ["mutation","position","wt_aa","mut_aa","chromosome","genomic_start","genomic_end","strand"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    mapped_muts = sum(1 for r in results if r["genomic_start"] != "NA")

    print("")
    print("=" * 65)
    print("  STEP 6 COMPLETE")
    print("=" * 65)
    print("  Mutations mapped:  " + str(mapped_muts) + " / " + str(len(results)))
    print("  Chromosome:        " + chrom)
    print("  Strand:            " + str(strand))

    for hs in ["R175H","R248W","R273H"]:
        for r in results:
            if r["mutation"] == hs:
                print("  " + hs + " -> chr" + str(r["chromosome"]) + ":" + str(r["genomic_start"]) + "-" + str(r["genomic_end"]))
                break

    print("")
    print("  Output: " + out_path)
    print("  Next: python code/step_07_gnomad_query.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
