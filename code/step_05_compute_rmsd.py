#!/usr/bin/env python3
import os, sys, csv, glob
from Bio.PDB import PDBParser, Superimposer
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(__file__))
from config import *

def extract_ca_atoms(structure):
    cas = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    cas.append(residue["CA"])
    return cas

def extract_plddt(pdb_path):
    plddts = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bfactor = float(line[60:66].strip())
                plddts.append(bfactor)
    return plddts

def main():
    print("=" * 65)
    print("  QuadSig Step 5: Compute RMSD + pLDDT Change")
    print("=" * 65)

    wt_pdb = os.path.join(DATA_STRUCT_WT, "wt_p53.pdb")
    if not os.path.exists(wt_pdb):
        print("[ERROR] WT PDB not found: " + wt_pdb)
        sys.exit(1)

    mut_path = os.path.join(DATA_INTER, "mutations.txt")
    with open(mut_path) as f:
        all_mutations = [line.strip() for line in f if line.strip()]
    print("[Step 5] Total mutations: " + str(len(all_mutations)))

    parser = PDBParser(QUIET=True)
    wt_struct = parser.get_structure("wt", wt_pdb)
    wt_cas = extract_ca_atoms(wt_struct)
    wt_plddt = extract_plddt(wt_pdb)
    print("[Step 5] WT CA atoms: " + str(len(wt_cas)))
    print("[Step 5] WT pLDDT values: " + str(len(wt_plddt)))
    print("[Step 5] WT avg pLDDT: " + str(round(sum(wt_plddt)/len(wt_plddt), 2)))

    results = []
    errors = []
    sup = Superimposer()

    print("[Step 5] Processing mutants...")
    for i, mutation in enumerate(all_mutations):
        mut_pdb = os.path.join(DATA_STRUCT_MUT, mutation + ".pdb")
        if not os.path.exists(mut_pdb):
            errors.append(mutation)
            continue

        try:
            mut_struct = parser.get_structure(mutation, mut_pdb)
            mut_cas = extract_ca_atoms(mut_struct)
            mut_plddt = extract_plddt(mut_pdb)

            pos = int(mutation[1:-1])
            pos_idx = pos - 1

            if len(mut_cas) == len(wt_cas):
                sup.set_atoms(wt_cas, mut_cas)
                rmsd = round(sup.rms, 4)
            else:
                rmsd = -1

            plddt_wt_pos = wt_plddt[pos_idx] if pos_idx < len(wt_plddt) else 0
            plddt_mut_pos = mut_plddt[pos_idx] if pos_idx < len(mut_plddt) else 0
            plddt_change = round(plddt_mut_pos - plddt_wt_pos, 4)

            wt_mean = sum(wt_plddt) / len(wt_plddt) if wt_plddt else 0
            mut_mean = sum(mut_plddt) / len(mut_plddt) if mut_plddt else 0

            start = max(0, pos_idx - 5)
            end = min(len(wt_plddt), pos_idx + 6)
            neigh_wt = wt_plddt[start:end]
            neigh_mut = mut_plddt[start:end] if len(mut_plddt) >= end else []
            neigh_wt_mean = round(sum(neigh_wt)/len(neigh_wt), 4) if neigh_wt else 0
            neigh_mut_mean = round(sum(neigh_mut)/len(neigh_mut), 4) if neigh_mut else 0

            results.append({
                "mutation": mutation,
                "position": pos,
                "wt_aa": mutation[0],
                "mut_aa": mutation[-1],
                "rmsd": rmsd,
                "plddt_wt": round(plddt_wt_pos, 4),
                "plddt_mut": round(plddt_mut_pos, 4),
                "plddt_change": plddt_change,
                "plddt_wt_mean": round(wt_mean, 4),
                "plddt_mut_mean": round(mut_mean, 4),
                "neighborhood_wt": neigh_wt_mean,
                "neighborhood_mut": neigh_mut_mean,
            })

        except Exception as e:
            errors.append(mutation + ": " + str(e))

        if (i + 1) % 500 == 0:
            print("    " + str(i + 1) + "/" + str(len(all_mutations)) + " processed")

    print("[Step 5] Processed: " + str(len(results)) + " | Errors: " + str(len(errors)))

    out_path = os.path.join(DATA_INTER, "structural_scores.csv")
    fields = ["mutation","position","wt_aa","mut_aa","rmsd","plddt_wt","plddt_mut","plddt_change","plddt_wt_mean","plddt_mut_mean","neighborhood_wt","neighborhood_mut"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    print("[Step 5] Saved: " + out_path)

    rmsds = [r["rmsd"] for r in results if r["rmsd"] >= 0]
    changes = [r["plddt_change"] for r in results]

    print("")
    print("=" * 65)
    print("  STEP 5 COMPLETE")
    print("=" * 65)
    print("  Mutations scored:  " + str(len(results)))
    print("  RMSD range:        " + str(round(min(rmsds), 3)) + " - " + str(round(max(rmsds), 3)))
    print("  RMSD mean:         " + str(round(sum(rmsds)/len(rmsds), 3)))
    print("  pLDDT change range: " + str(round(min(changes), 3)) + " to " + str(round(max(changes), 3)))
    print("  Errors:            " + str(len(errors)))

    hotspots = ["R175H","G245S","R248W","R248Q","R249S","R273H","R273C"]
    print("  Cancer hotspot RMSD:")
    for hs in hotspots:
        for r in results:
            if r["mutation"] == hs:
                print("    " + hs + ": RMSD=" + str(r["rmsd"]) + " pLDDT_change=" + str(r["plddt_change"]))
                break

    print("")
    print("  Output: " + out_path + " (" + str(os.path.getsize(out_path)) + " bytes)")
    print("  Next: python code/step_06_genomic_coords.py")
    print("=" * 65)

if __name__ == "__main__":
    main()
