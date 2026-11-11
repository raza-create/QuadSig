#!/usr/bin/env python3
"""Recompute per-position conservation from ONE orthologue per vertebrate species.

Input : data/raw/p53_msa.fasta (Clustal Omega alignment of the human sequence and every
        Ensembl Compara TP53 orthologue; 254 rows, several species with more than one row)
Rule  : for each non-human species keep the row with the highest identity to human p53 over
        the 393 human columns (gaps count as mismatches; ties -> longer sequence);
        drop non-vertebrate species (Ensembl's tunicates/invertebrates/yeast).
Score : conservation(i) = fraction of the retained orthologues (human excluded) carrying the
        human residue in the column of human position i, gaps excluded; rounded to 4 dp.
Output: data/intermediate/conservation_scores.csv            (same columns as before)
        data/intermediate/ortholog_selection.csv             (every alignment row, kept or not, and why)
        data/raw/p53_msa_one_per_species.fasta               (the rows used: human + orthologues)
The previous conservation file is kept as data/intermediate/conservation_scores_msa254.csv.
Usage (from ~/QuadSig):  python3 recompute_conservation.py
"""
import os, re, shutil, sys
import numpy as np
import pandas as pd

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
MSA = os.path.join(ROOT, "data/raw/p53_msa.fasta")
CONS = os.path.join(ROOT, "data/intermediate/conservation_scores.csv")
OLD = os.path.join(ROOT, "data/intermediate/conservation_scores_msa254.csv")
P04637 = ("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ"
          "KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRV"
          "EGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEP"
          "HHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
NONVERT = ("ciona_", "drosophila_", "caenorhabditis_", "saccharomyces_")   # non-vertebrates in Ensembl
DOMAINS = {"TAD": (1, 67), "PRR": (68, 98), "DBD": (102, 292), "TET": (323, 356), "CTD": (364, 393)}  # = config.py
ZINC = (176, 179, 238, 242)


def fasta(path):
    out, name, buf = [], None, []
    for line in open(path, errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif name is not None:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


msa = [(n.split()[0], s.upper()) for n, s in fasta(MSA)]
hum = [i for i, (n, s) in enumerate(msa) if s.replace("-", "").replace(".", "") == P04637]
hi = next((i for i in hum if msa[i][0].lower().startswith("human")), hum[0] if hum else None)
if hi is None:
    sys.exit("STOPPED: no alignment row matches human P04637 exactly (nothing was changed)")
H = msa[hi][1]
cols = [c for c, a in enumerate(H) if a not in "-."]
assert len(cols) == 393
dbd = cols[101:292]

rows = []
for i, (n, s) in enumerate(msa):
    sp = "homo_sapiens" if i == hi else re.sub(r"_P53(_\d+)?$", "", n, flags=re.I).lower()
    rows.append(dict(row=i, name=n, species=sp, length=len(s.replace("-", "").replace(".", "")),
                     identity=round(sum(s[c] == H[c] for c in cols) / 393, 4),
                     dbd_identity=round(sum(s[c] == H[c] for c in dbd) / len(dbd), 4)))
sel = pd.DataFrame(rows)
sel["n_rows_species"] = sel.groupby("species").row.transform("size")
sel["vertebrate"] = ~sel.species.str.startswith(NONVERT)
sel["kept"] = False
sel["reason"] = ""
sel.loc[sel.row == hi, "reason"] = "human reference"
cand = sel[(sel.row != hi) & (sel.species != "homo_sapiens") & sel.vertebrate]
best = cand.sort_values(["species", "identity", "length"], ascending=[True, False, False]).groupby("species").head(1)
sel.loc[best.index, "kept"] = True
sel.loc[best.index, "reason"] = "highest identity to human for this species"
sel.loc[(sel.row != hi) & ~sel.vertebrate, "reason"] = "non-vertebrate species"
sel.loc[(sel.reason == "") & (sel.species == "homo_sapiens"), "reason"] = "extra human row"
sel.loc[sel.reason == "", "reason"] = "another row of this species is more similar to human"
kept = sel[sel.kept]

# ---- conservation from the retained orthologues (human excluded), gaps excluded ----
K = [msa[i][1] for i in kept.row]
cons, nres = [], []
for c in cols:
    res = [s[c] for s in K if s[c] not in "-."]
    nres.append(len(res))
    cons.append(round(sum(r == H[c] for r in res) / len(res), 4) if res else np.nan)

# ---- write (old file kept once) ----
if not os.path.exists(OLD):
    shutil.copy2(CONS, OLD)
old = pd.read_csv(OLD)
new = old.copy()
new["conservation_score"] = new.position.map(dict(zip(range(1, 394), cons)))
new.to_csv(CONS, index=False)
sel.drop(columns="row").to_csv(os.path.join(ROOT, "data/intermediate/ortholog_selection.csv"), index=False)
with open(os.path.join(ROOT, "data/raw/p53_msa_one_per_species.fasta"), "w") as f:
    for i in [hi] + list(kept.row):
        f.write(f">{msa[i][0]}\n{msa[i][1]}\n")

# ---- report ----
print(f"alignment rows: {len(msa)}; species: {sel.species.nunique()} (incl. human)")
print(f"removed non-vertebrate species: {sorted(sel.loc[~sel.vertebrate, 'species'].unique()) or 'none'}")
print(f"ORTHOLOGUES USED: {len(kept)} (one per vertebrate species) + human = {len(kept) + 1} sequences")
print(f"  rows not used: {len(sel) - len(kept) - 1}  "
      f"(extra rows of the same species: {(sel.reason == 'another row of this species is more similar to human').sum()}, "
      f"non-vertebrate: {(~sel.vertebrate).sum()})")
print(f"  identity of used orthologues to human: min {kept.identity.min():.2f}, median {kept.identity.median():.2f}; "
      f"lowest: " + ", ".join(f"{r.species} {r.identity:.2f}" for r in kept.nsmallest(4, 'identity').itertuples()))
extra = [c for c in old.columns if c not in ("position", "conservation_score")]
if extra:
    print(f"  NOTE: conservation_scores.csv also has columns {extra}; they were NOT recomputed")
print(f"  positions with no orthologue residue (score missing): {int(np.isnan(cons).sum())}; "
      f"orthologues per position: min {min(nres)}, median {int(np.median(nres))}")
o, n_ = old.set_index("position").conservation_score, pd.Series(cons, index=range(1, 394))
print(f"\nconservation  old (all {len(msa)} rows)  ->  new (one per species)")
print(f"  range        {o.min():.3f}-{o.max():.3f}  ->  {np.nanmin(cons):.3f}-{np.nanmax(cons):.3f}")
print(f"  mean |change| {np.nanmean(np.abs(o.values - n_.values)):.3f}; max |change| {np.nanmax(np.abs(o.values - n_.values)):.3f}; "
      f"Spearman old vs new {pd.Series(o.values).corr(pd.Series(n_.values), method='spearman'):.3f}")
for d, (a, b) in DOMAINS.items():
    print(f"  {d:4s} mean   {o.loc[a:b].mean():.2f}  ->  {n_.loc[a:b].mean():.2f}")
print("  zinc residues " + ", ".join(f"{P04637[p-1]}{p} {o.loc[p]:.2f}->{n_.loc[p]:.2f}" for p in ZINC))
print(f"\nwrote {CONS}\n      data/intermediate/ortholog_selection.csv, data/raw/p53_msa_one_per_species.fasta"
      f"\nkept  {OLD} (previous scores)")
