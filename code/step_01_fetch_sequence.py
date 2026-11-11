#!/usr/bin/env python3
import os, sys, json, csv, requests

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def write_fasta(filepath, header, sequence):
    with open(filepath, "w") as f:
        f.write(">" + header + "\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

def main():
    print("=" * 65)
    print("  QuadSig Step 1: Fetch p53 Wild-Type Sequence")
    print("=" * 65)

    # Download FASTA from UniProt
    url = UNIPROT_API + "/" + UNIPROT_ID + ".fasta"
    print("[Step 1] Downloading from: " + url)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    lines = resp.text.strip().split("\n")
    seq = "".join(l for l in lines if not l.startswith(">"))
    print("[Step 1] Downloaded sequence: " + str(len(seq)) + " aa")

    # Download full JSON record for annotations
    url2 = UNIPROT_API + "/" + UNIPROT_ID + ".json"
    print("[Step 1] Downloading annotations: " + url2)
    resp2 = requests.get(url2, timeout=60)
    resp2.raise_for_status()
    record = resp2.json()
    json_path = os.path.join(DATA_RAW, "p53_uniprot.json")
    with open(json_path, "w") as f:
        json.dump(record, f, indent=2)
    print("[Step 1] Saved UniProt JSON: " + json_path)

    # Validate length
    if len(seq) != SEQUENCE_LEN:
        print("[ERROR] Expected " + str(SEQUENCE_LEN) + " aa, got " + str(len(seq)))
        sys.exit(1)
    print("[Step 1] Sequence length verified: " + str(len(seq)) + " aa")

    # Verify zinc residues
    print("[Step 1] Checking zinc-coordination residues:")
    for pos in sorted(ZINC_RESIDUES.keys()):
        expected = ZINC_RESIDUES[pos]
        actual = seq[pos - 1]
        if actual == expected:
            print("    " + expected + str(pos) + ": PASS")
        else:
            print("    " + expected + str(pos) + ": FAIL (found " + actual + ")")
            sys.exit(1)
    print("[Step 1] All zinc residues verified")

    # Save FASTA files
    write_fasta(os.path.join(DATA_RAW, "p53_sequence.fasta"), "sp|" + UNIPROT_ID + "|P53_HUMAN Cellular tumor antigen p53", seq)
    write_fasta(os.path.join(DATA_FASTA, "wt_p53.fasta"), "WT_P53", seq)
    print("[Step 1] Saved FASTA files")

    # Save plain text sequence
    txt_path = os.path.join(DATA_RAW, "p53_sequence.txt")
    with open(txt_path, "w") as f:
        f.write(seq)
    print("[Step 1] Saved plain sequence: " + txt_path)

    # Extract annotations from UniProt JSON
    annotations = []
    for feat in record.get("features", []):
        loc = feat.get("location", {})
        s = loc.get("start", {}).get("value")
        e = loc.get("end", {}).get("value")
        if s is None or e is None:
            continue
        for p in range(s, e + 1):
            if 1 <= p <= len(seq):
                annotations.append({
                    "position": p,
                    "wt_aa": seq[p - 1],
                    "feature_type": feat.get("type", ""),
                    "description": feat.get("description", "")
                })

    ann_path = os.path.join(DATA_RAW, "p53_annotations.csv")
    with open(ann_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["position", "wt_aa", "feature_type", "description"])
        w.writeheader()
        w.writerows(annotations)
    print("[Step 1] Saved " + str(len(annotations)) + " annotations: " + ann_path)

    # Summary
    print("")
    print("=" * 65)
    print("  STEP 1 COMPLETE")
    print("=" * 65)
    print("  Sequence length:  " + str(len(seq)) + " aa")
    print("  First 20 aa:      " + seq[:20] + "...")
    print("  Last 20 aa:       ..." + seq[-20:])
    print("  Zinc residues:    C176 H179 C238 C242 — all verified")
    print("  Annotations:      " + str(len(annotations)) + " entries from UniProt")
    print("")
    print("  Files created:")
    created = [
        "data/raw/p53_uniprot.json",
        "data/raw/p53_sequence.fasta",
        "data/raw/p53_sequence.txt",
        "data/raw/p53_annotations.csv",
        "data/fasta/wt_p53.fasta"
    ]
    for p in created:
        full = os.path.join(PROJECT_ROOT, p)
        if os.path.exists(full):
            sz = os.path.getsize(full)
            print("    [OK] " + p + "  (" + str(sz) + " bytes)")
        else:
            print("    [MISSING] " + p)
    print("")
    print("  Next: python code/step_02_generate_mutants.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
