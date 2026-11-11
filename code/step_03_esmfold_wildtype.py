#!/usr/bin/env python3
import os, sys, time, gc, torch

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 65)
    print("  QuadSig Step 3: ESMFold on Wild-Type p53")
    print("=" * 65)

    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_name(0)
        vram = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        print("[Step 3] GPU: " + gpu_name + " (" + str(round(vram, 1)) + " GB)")

    seq_path = os.path.join(DATA_RAW, "p53_sequence.txt")
    with open(seq_path) as f:
        seq = f.read().strip()
    print("[Step 3] Loaded WT sequence: " + str(len(seq)) + " aa")

    from transformers import EsmForProteinFolding, AutoTokenizer
    from transformers.models.esm.openfold_utils.protein import to_pdb, Protein
    import numpy as np

    print("[Step 3] Loading ESMFold model...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    gc.collect()
    model.eval()
    model.esm = model.esm.half()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        model = model.cuda()
    model.trunk.set_chunk_size(64)
    print("[Step 3] Model loaded in " + str(round(time.time() - t0, 1)) + " sec")

    print("[Step 3] Predicting WT p53 structure...")
    t1 = time.time()
    tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
    if torch.cuda.is_available():
        tokenized = {k: v.cuda() for k, v in tokenized.items()}
    with torch.no_grad():
        output = model(**tokenized)
    pred_time = round(time.time() - t1, 1)
    print("[Step 3] Prediction done in " + str(pred_time) + " sec")

    # Extract pLDDT - handle any shape
    raw_plddt = output["plddt"][0].cpu().float()
    if raw_plddt.dim() == 0:
        plddt_values = [raw_plddt.item()]
    elif raw_plddt.dim() == 1:
        plddt_values = raw_plddt.tolist()
    else:
        plddt_values = raw_plddt.mean(dim=-1).tolist()
    plddt_np = np.array(plddt_values)

    # Build PDB
    pdb_string = to_pdb(Protein(
        aatype=tokenized["input_ids"][0].cpu().numpy(),
        atom_positions=output["positions"][-1][0].cpu().float().numpy(),
        atom_mask=output["atom37_atom_exists"][0].cpu().float().numpy(),
        residue_index=np.arange(len(seq), dtype=np.int32),
        b_factors=np.tile(plddt_np[:, None], (1, 37)),
    ))
    pdb_path = os.path.join(DATA_STRUCT_WT, "wt_p53.pdb")
    with open(pdb_path, "w") as f:
        f.write(pdb_string)
    print("[Step 3] Saved PDB: " + pdb_path)

    # Stats
    n = min(len(plddt_values), len(seq))
    avg_p = sum(plddt_values[:n]) / n
    min_p = min(plddt_values[:n])
    max_p = max(plddt_values[:n])
    print("[Step 3] pLDDT: " + str(n) + " residues")
    print("[Step 3] Average pLDDT: " + str(round(avg_p, 2)))
    print("[Step 3] Min pLDDT: " + str(round(min_p, 2)))
    print("[Step 3] Max pLDDT: " + str(round(max_p, 2)))

    # Save pLDDT CSV
    plddt_path = os.path.join(DATA_INTER, "wt_plddt.csv")
    with open(plddt_path, "w") as f:
        f.write("position,wt_aa,plddt\n")
        for i in range(n):
            f.write(str(i+1) + "," + seq[i] + "," + str(round(plddt_values[i], 2)) + "\n")
    print("[Step 3] Saved pLDDT: " + plddt_path)

    print("[Step 3] Zinc residue pLDDT:")
    for pos in sorted(ZINC_RESIDUES.keys()):
        if pos <= n:
            val = plddt_values[pos - 1]
            status = "HIGH" if val > 70 else "LOW"
            print("    " + ZINC_RESIDUES[pos] + str(pos) + ": " + str(round(val, 2)) + " (" + status + ")")

    print("")
    print("=" * 65)
    print("  STEP 3 COMPLETE")
    print("=" * 65)
    print("  Prediction time:  " + str(pred_time) + " sec")
    print("  PDB file:         " + pdb_path + " (" + str(os.path.getsize(pdb_path)) + " bytes)")
    print("  Residues:         " + str(n))
    print("  Avg pLDDT:        " + str(round(avg_p, 2)))
    print("")
    print("  Files created:")
    print("    [OK] " + pdb_path)
    print("    [OK] " + plddt_path)
    print("")
    print("  Next: python code/step_04_esmfold_mutants.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
