#!/usr/bin/env python3
import os, sys, time, glob, torch, numpy as np
from transformers import EsmForProteinFolding, AutoTokenizer
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def output_to_pdb(output):
    final_pos = output["positions"][-1].cpu().numpy()
    output_np = {k: v.cpu().numpy() for k, v in output.items() if k != "positions"}
    mask = output_np["atom37_atom_exists"]
    pdb_str = to_pdb(Protein(
        aatype=output_np["aatype"][0],
        atom_positions=final_pos[0],
        atom_mask=mask[0],
        residue_index=output_np["residue_index"][0] + 1,
        b_factors=output_np["plddt"][0],
    ))
    return pdb_str

def main():
    print("=" * 65)
    print("  QuadSig Step 4: ESMFold on All 7,467 Mutant Sequences")
    print("=" * 65)

    gpu = torch.cuda.get_device_name(0)
    vram = torch.cuda.get_device_properties(0).total_memory / (1024**3)
    print("[Step 4] GPU: " + gpu + " (" + str(round(vram, 1)) + " GB)")

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 4] Total mutations: " + str(len(all_mutations)))

    done = set()
    for pf in glob.glob(os.path.join(DATA_STRUCT_MUT, "*.pdb")):
        done.add(os.path.basename(pf).replace(".pdb", ""))
    remaining = [m for m in all_mutations if m not in done]
    print("[Step 4] Already completed: " + str(len(done)))
    print("[Step 4] Remaining: " + str(len(remaining)))

    if len(remaining) == 0:
        print("[Step 4] All done! Next: python code/step_05_compute_rmsd.py")
        return

    seq_path = os.path.join(DATA_RAW, "p53_sequence.txt")
    with open(seq_path) as f:
        wt_seq = f.read().strip()

    print("[Step 4] Loading ESMFold model...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model.esm = model.esm.half()
    model.trunk.set_chunk_size(64)
    model = model.eval().cuda()
    torch.cuda.empty_cache()
    print("[Step 4] Model loaded in " + str(round(time.time() - t0, 1)) + " sec")
    print("[Step 4] Free VRAM: " + str(round(torch.cuda.mem_get_info()[0]/1024**3, 1)) + " GB")

    log_path = os.path.join(LOGS_DIR, "step04_progress.log")
    total_done = len(done)
    total_all = len(all_mutations)
    start_time = time.time()

    print("[Step 4] Starting predictions...")
    print("[Step 4] Safe to Ctrl+C. Re-run to resume.")
    print("")

    for i, mutation in enumerate(remaining):
        mut_aa = mutation[-1]
        pos = int(mutation[1:-1])
        pos_idx = pos - 1
        mut_seq = wt_seq[:pos_idx] + mut_aa + wt_seq[pos_idx + 1:]

        t1 = time.time()
        try:
            inputs = tokenizer([mut_seq], return_tensors="pt", add_special_tokens=False)
            inputs = {k: v.cuda() for k, v in inputs.items()}

            with torch.no_grad():
                output = model(**inputs)

            pdb_str = output_to_pdb(output)
            pdb_path = os.path.join(DATA_STRUCT_MUT, mutation + ".pdb")
            with open(pdb_path, "w") as f:
                f.write(pdb_str)

            del output, inputs
            pred_time = round(time.time() - t1, 1)
            total_done += 1

            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta_sec = (len(remaining) - i - 1) / rate if rate > 0 else 0
            eta_hrs = round(eta_sec / 3600, 1)
            eta_days = round(eta_sec / 86400, 1)

            print("    [" + str(total_done) + "/" + str(total_all) + "] " + mutation + " " + str(pred_time) + "s | ETA: " + str(eta_hrs) + "h (" + str(eta_days) + "d)")

            with open(log_path, "a") as f:
                f.write(mutation + "," + str(pred_time) + "," + str(total_done) + "/" + str(total_all) + "\n")

        except Exception as e:
            print("    [ERROR] " + mutation + ": " + str(e))
            with open(log_path, "a") as f:
                f.write(mutation + ",ERROR," + str(e) + "\n")
            torch.cuda.empty_cache()
            continue

        if (i + 1) % 25 == 0:
            torch.cuda.empty_cache()

    total_pdb = len(glob.glob(os.path.join(DATA_STRUCT_MUT, "*.pdb")))
    elapsed_total = round((time.time() - start_time) / 3600, 1)

    print("")
    print("=" * 65)
    print("  STEP 4 COMPLETE")
    print("=" * 65)
    print("  PDB files: " + str(total_pdb) + " / " + str(total_all))
    print("  Time: " + str(elapsed_total) + " hours")
    if total_pdb == total_all:
        print("  ALL 7,467 DONE. Next: python code/step_05_compute_rmsd.py")
    else:
        print("  Re-run to continue from checkpoint")
    print("=" * 65)

if __name__ == "__main__":
    main()
