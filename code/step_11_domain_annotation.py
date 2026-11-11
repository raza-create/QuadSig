#!/usr/bin/env python3
import os, sys, csv

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 11: Functional Domain Annotation")
    print("=" * 65)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 11] Total mutations: " + str(len(all_mutations)))

    def get_domain(pos):
        for code, (start, end, name) in DOMAINS.items():
            if start <= pos <= end:
                return code, name
        return "LINKER", "Linker region"

    results = []
    domain_counts = {}
    for mutation in all_mutations:
        pos = int(mutation[1:-1])
        domain_code, domain_name = get_domain(pos)
        results.append({
            "mutation": mutation,
            "position": pos,
            "wt_aa": mutation[0],
            "mut_aa": mutation[-1],
            "domain_code": domain_code,
            "domain_name": domain_name,
        })
        domain_counts[domain_code] = domain_counts.get(domain_code, 0) + 1

    out_path = os.path.join(DATA_INTER, "domain_annotations.csv")
    fields = ["mutation","position","wt_aa","mut_aa","domain_code","domain_name"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    print("[Step 11] Domain distribution:")
    for code, (start, end, name) in DOMAINS.items():
        count = domain_counts.get(code, 0)
        positions = end - start + 1
        print("    " + code + " (" + str(start) + "-" + str(end) + "): " + str(positions) + " positions, " + str(count) + " mutations — " + name)
    linker = domain_counts.get("LINKER", 0)
    print("    LINKER: " + str(linker) + " mutations")

    # Hotspot domains
    print("[Step 11] Cancer hotspot domains:")
    for hs in ["R175H","G245S","R248W","R248Q","R249S","R273H","R273C"]:
        for r in results:
            if r["mutation"] == hs:
                print("    " + hs + " -> " + r["domain_code"] + " (" + r["domain_name"] + ")")
                break

    print("")
    print("=" * 65)
    print("  STEP 11 COMPLETE")
    print("=" * 65)
    print("  Mutations annotated: " + str(len(results)))
    print("  Output: " + out_path)
    print("  Next: python code/step_12_zinc_flagging.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
