#!/usr/bin/env python3
import os, sys, requests, time

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 8: Download p53 MSA Across Species")
    print("=" * 65)

    # Use UniProt cluster/homolog search to get ortholog sequences
    # Then align them using EBI Clustal Omega API

    # Step 1: Get p53 orthologs from Ensembl
    print("[Step 8] Fetching p53 orthologs from Ensembl...")
    url = ENSEMBL_API + "/homology/id/human/ENSG00000141510?type=orthologues;sequence=protein;format=full"
    resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=120)
    resp.raise_for_status()
    data = resp.json()

    homologies = data["data"][0]["homologies"]
    print("[Step 8] Total orthologs found: " + str(len(homologies)))

    # Collect sequences
    sequences = {}
    # Add human WT first
    seq_path = os.path.join(DATA_RAW, "p53_sequence.txt")
    with open(seq_path) as f:
        human_seq = f.read().strip()
    sequences["Human_P53"] = human_seq

    for h in homologies:
        species = h["target"]["species"].replace(" ", "_")
        target_seq = h["target"].get("align_seq", "")
        if not target_seq:
            target_seq = h["target"].get("seq", "")
        if target_seq and len(target_seq) > 100:
            # Remove gaps for raw sequence
            clean_seq = target_seq.replace("-", "")
            if len(clean_seq) > 100:
                name = species + "_P53"
                sequences[name] = clean_seq

    print("[Step 8] Sequences with length > 100 aa: " + str(len(sequences)))

    # Save raw orthologs FASTA
    ortho_path = os.path.join(DATA_RAW, "p53_orthologs.fasta")
    with open(ortho_path, "w") as f:
        for name, seq in sequences.items():
            f.write(">" + name + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    print("[Step 8] Saved orthologs: " + ortho_path)

    # Step 2: Submit to Clustal Omega for MSA
    print("[Step 8] Submitting to EBI Clustal Omega for alignment...")
    print("[Step 8] (This may take 1-5 minutes)")

    fasta_str = ""
    for name, seq in sequences.items():
        fasta_str += ">" + name + "\n" + seq + "\n"

    submit_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    params = {
        "email": "quadsig@example.com",
        "sequence": fasta_str,
        "stype": "protein",
        "outfmt": "fa",
    }
    resp = requests.post(submit_url, data=params, timeout=120)
    if resp.status_code != 200:
        print("[Step 8] Clustal Omega submit failed: " + str(resp.status_code))
        print("[Step 8] Saving unaligned sequences as fallback...")
        msa_path = os.path.join(DATA_RAW, "p53_msa.fasta")
        with open(msa_path, "w") as f:
            f.write(fasta_str)
        print("[Step 8] Saved unaligned FASTA: " + msa_path)
    else:
        job_id = resp.text.strip()
        print("[Step 8] Job submitted: " + job_id)

        # Poll for completion
        status_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/" + job_id
        for attempt in range(60):
            time.sleep(10)
            resp = requests.get(status_url, timeout=30)
            status = resp.text.strip()
            if status == "FINISHED":
                print("[Step 8] Alignment complete!")
                break
            elif status == "FAILURE" or status == "ERROR":
                print("[Step 8] Alignment failed: " + status)
                break
            else:
                if attempt % 3 == 0:
                    print("[Step 8] Status: " + status + " (waiting...)")
        
        if status == "FINISHED":
            result_url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/" + job_id + "/aln-fasta"
            resp = requests.get(result_url, timeout=120)
            msa_path = os.path.join(DATA_RAW, "p53_msa.fasta")
            with open(msa_path, "w") as f:
                f.write(resp.text)
            print("[Step 8] Saved MSA: " + msa_path)
        else:
            print("[Step 8] Using unaligned sequences as fallback")
            msa_path = os.path.join(DATA_RAW, "p53_msa.fasta")
            with open(msa_path, "w") as f:
                f.write(fasta_str)

    # Count sequences in MSA
    msa_path = os.path.join(DATA_RAW, "p53_msa.fasta")
    seq_count = 0
    with open(msa_path) as f:
        for line in f:
            if line.startswith(">"):
                seq_count += 1

    species_examples = []
    with open(msa_path) as f:
        for line in f:
            if line.startswith(">") and len(species_examples) < 10:
                species_examples.append(line.strip()[1:])

    print("")
    print("=" * 65)
    print("  STEP 8 COMPLETE")
    print("=" * 65)
    print("  Sequences in MSA:  " + str(seq_count))
    print("  File:              " + msa_path + " (" + str(os.path.getsize(msa_path)) + " bytes)")
    print("  Example species:")
    for s in species_examples:
        print("    " + s)
    print("")
    print("  Next: python code/step_09_conservation_scores.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
