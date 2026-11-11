#!/usr/bin/env python3
"""Supplementary Data workbook for the QuadSig v2 manuscript.

Builds results_v2/supplementary/QuadSig_Supplementary_Data.xlsx from files that already exist:
  results_v2/  results.csv, parameters.json, cv_folds.csv, clinvar_records.csv, summary.txt
  data/intermediate/  mutations.txt and the per-signal CSV files
  data/raw/p53_msa.fasta + data/intermediate/ortholog_selection.csv (recompute_conservation.py)
Nothing is recomputed or changed; every sheet is a copy of a pipeline output. Sheet 2 is
"mutations" and sheet 4 is "Conservation_Scores", as cited in the manuscript. Author/creator
metadata is left blank for double-anonymous review.
Usage (from ~/QuadSig):  python3 make_supplementary.py          [or: python3 make_supplementary.py ROOT]
"""
import json, math, os, re, sys
import pandas as pd

try:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
except ImportError:
    sys.exit("openpyxl is not installed. Install it safely with:  sudo apt install python3-openpyxl\n"
             "(nothing was written)")

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
RES, INTER = os.path.join(ROOT, "results_v2"), os.path.join(ROOT, "data", "intermediate")
OUTDIR = os.path.join(RES, "supplementary")
OUT = os.path.join(OUTDIR, "QuadSig_Supplementary_Data.xlsx")
problems = []


def need(path):
    if not os.path.exists(path):
        sys.exit(f"MISSING: {path}  (nothing was written)")
    return path


def rcsv(path):
    return pd.read_csv(need(path), dtype=str, keep_default_na=False)


# ---------------- load (as text, so values are copied exactly as the pipeline wrote them) ----
res = rcsv(os.path.join(RES, "results.csv"))
params = json.load(open(need(os.path.join(RES, "parameters.json"))))
folds = rcsv(os.path.join(RES, "cv_folds.csv"))
records = rcsv(os.path.join(RES, "clinvar_records.csv"))
summary = open(need(os.path.join(RES, "summary.txt"))).read()
muts = [l.strip() for l in open(need(os.path.join(INTER, "mutations.txt"))) if l.strip()]
struct = rcsv(os.path.join(INTER, "structural_scores.csv"))
pop = rcsv(os.path.join(INTER, "population_scores.csv"))
cons = rcsv(os.path.join(INTER, "conservation_scores.csv"))
phys = rcsv(os.path.join(INTER, "physicochemical_scores.csv"))
dom = rcsv(os.path.join(INTER, "domain_annotations.csv"))
zn = rcsv(os.path.join(INTER, "zinc_flags.csv"))
gof = rcsv(os.path.join(INTER, "gof_flags.csv"))
m = None
for f_ in (os.path.join(RES, "summary.txt"), os.path.join(RES, "clinvar_retrieval_date.txt"), os.path.join(ROOT, "results_v2_msa254", "summary.txt")):
    if os.path.exists(f_):          # an --offline rerun reuses the records; the fetch date is in the first run
        m = m or re.search(r"fetched (\d{4}-\d{2}-\d{2})", open(f_).read())
clinvar_date = m.group(1) if m else "DATE NOT FOUND"
if not m:
    problems.append("ClinVar retrieval date not found in summary.txt; fix it in the README sheet")


# ---------------- orthologue alignment: all 254 rows, with the one-per-species selection ----------
def parse_fasta(path):
    out, name, buf = [], None, []
    for line in open(need(path), errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip().split()[0], []
        elif name is not None:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


msa_rows = parse_fasta(os.path.join(ROOT, "data", "raw", "p53_msa.fasta"))
selection = rcsv(os.path.join(INTER, "ortholog_selection.csv"))
if [n for n, _ in msa_rows] != list(selection.name):
    problems.append("ortholog_selection.csv rows do not match p53_msa.fasta (rerun recompute_conservation.py)")
orth = selection.copy()
orth["aligned_sequence"] = [s_ for _, s_ in msa_rows][:len(orth)]
n_kept = int((orth.kept == "True").sum())
print(f"Orthologue alignment: {len(msa_rows)} rows; used for conservation: {n_kept} + human")

# ---------------- consistency checks (report only; nothing is changed) -----------------------
if len(res) != len(muts) or set(res.mutation) != set(muts):
    problems.append(f"results.csv ({len(res)}) and mutations.txt ({len(muts)}) do not list the same mutations")
for name, df in [("structural_scores", struct), ("population_scores", pop), ("physicochemical_scores", phys),
                 ("domain_annotations", dom), ("zinc_flags", zn), ("gof_flags", gof)]:
    if len(df) != len(muts):
        problems.append(f"{name}.csv has {len(df)} rows (expected {len(muts)})")
if len(cons) != 393:
    problems.append(f"conservation_scores.csv has {len(cons)} rows (expected 393)")
lab = res[res.clinvar_label.isin(["pathogenic", "benign"])]

# ---------------- sheets -----------------------------------------------------------------------
mut_sheet = res[["mutation", "position", "wt_aa", "mut_aa", "domain_code", "domain_name"]].copy()
annot = (dom.merge(zn, on="mutation", how="outer").merge(gof, on="mutation", how="outer"))
clin = res.loc[res.clinvar_label != "", ["mutation", "clinvar_label", "n_records", "clinvar_significance",
                                         "quadsig_score", "category", "clinvar_agreement"]]
param_rows = [("weights." + k, v) for k, v in params["weights"].items()]
param_rows += [("domain_term." + k, v) for k, v in params["domain_add"].items()]
param_rows += [("threshold_low (lowest fold threshold)", params["threshold_low"]),
               ("threshold_high (highest fold threshold)", params["threshold_high"]),
               ("youden_threshold_all_labelled (in-sample)", params["youden_all"]),
               ("fold_thresholds", ", ".join(f"{t:.4f}" for t in params["fold_thresholds"])),
               ("roc_auc", params["auc"])]
param_rows += [("cross_validated." + k, v) for k, v in params["cv"].items()]
param_rows += [("in_sample." + k, v) for k, v in params["in_sample"].items()]
param_df = pd.DataFrame(param_rows, columns=["parameter", "value"])

RESULT_COLS = {
    "mutation": "Variant (wild-type residue, position, mutant residue)",
    "position": "Residue position in UniProt P04637 (1-393)",
    "wt_aa": "Wild-type residue", "mut_aa": "Mutant residue",
    "quadsig_score": "QuadSig score (0-1), formula in Methods Section 2.3",
    "category": "likely_pathogenic / ambiguous / likely_benign (cross-validated thresholds)",
    "population_score": "Population (rarity) score from gnomAD v4 allele frequency (Table 1)",
    "allele_frequency": "Maximum gnomAD v4 allele frequency over exome and genome call sets (0 = absent)",
    "conservation_score": "Fraction of the vertebrate orthologues (one per species) carrying the human residue in that alignment column, gaps excluded",
    "rmsd": "C-alpha RMSD (Angstrom) of the ESMFold mutant versus the wild type",
    "rmsd_norm": "Logistic RMSD score, 1/(1+exp(-1.5(RMSD-2)))",
    "plddt_wt": "Wild-type pLDDT (0-100), from structural_scores.csv",
    "plddt_mut": "Mutant pLDDT (0-100), from structural_scores.csv",
    "plddt_change": "Mutant minus wild-type pLDDT",
    "plddt_penalty": "min(|pLDDT change|/30, 1) when pLDDT decreases, otherwise 0",
    "physicochemical_score": "RMS of normalised changes in charge, mass and hydrophobicity",
    "blosum62": "BLOSUM62 score of the substitution (used only for consistency check 6)",
    "domain_code": "TAD, PRR, DBD, TET, CTD or LINKER", "domain_name": "Domain name",
    "is_zinc_coord": "Substitution at C176, H179, C238 or C242 (annotation only; not scored)",
    "is_gof": "One of the seven curated GOF hotspots (annotation only; not scored)",
    "mechanism": "Signals supporting the score (text label)",
    "clinvar_label": "ClinVar protein-level consensus label (pathogenic / benign / uncertain; blank = not in ClinVar)",
    "n_records": "Number of ClinVar records for this protein change",
    "clinvar_significance": "ClinVar germline classifications of those records",
    "clinvar_agreement": "AGREE / DISAGREE for labelled variants (likely pathogenic vs pathogenic, likely benign vs benign; ambiguous counts as DISAGREE)",
}

sheets = [
    ("README", None),
    ("mutations", mut_sheet),
    ("Structural_Scores", struct),
    ("Conservation_Scores", cons),
    ("Population_Scores", pop),
    ("Physicochemical_Scores", phys),
    ("Annotations", annot),
    ("QuadSig_Results", res),
    ("ClinVar_Protein_Labels", clin),
    ("ClinVar_Records", records),
    ("CV_Folds", folds),
    ("Parameters", param_df),
    ("Orthologue_Alignment", orth),
]
DESC = {
    "mutations": "All 7,467 single amino-acid substitutions (393 positions x 19) with their domain.",
    "Structural_Scores": "Per-variant ESMFold metrics as written by the structural step (data/intermediate/structural_scores.csv).",
    "Conservation_Scores": "Per-position conservation score (393 positions), one orthologue per vertebrate species (see Orthologue_Alignment).",
    "Population_Scores": "Per-variant gnomAD v4 allele frequency and population score.",
    "Physicochemical_Scores": "Per-variant physicochemical score, as written by the pipeline (data/intermediate/physicochemical_scores.csv).",
    "Annotations": "Domain annotation, zinc-coordination flag and GOF hotspot flag with reported mechanism.",
    "QuadSig_Results": "Final table: QuadSig score, class, all signals and ClinVar label for every variant (column definitions below).",
    "ClinVar_Protein_Labels": "Every variant with a ClinVar germline classification: consensus label, records and QuadSig result.",
    "ClinVar_Records": f"All TP53 single-nucleotide variant records retrieved from ClinVar (E-utilities) on {clinvar_date}.",
    "CV_Folds": "Threshold and performance in each of the five stratified cross-validation folds (seed 42).",
    "Parameters": "Signal weights, domain terms, thresholds, ROC-AUC, cross-validated and in-sample performance.",
    "Orthologue_Alignment": "All rows of the Clustal Omega alignment (human + every Ensembl Compara TP53 orthologue sequence). "
                            "'kept' marks the one sequence per vertebrate species used for conservation (highest identity to human); "
                            "'reason' says why each other row was not used.",
}

# ---------------- write ------------------------------------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)
F, FB = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
FT = Font(name="Arial", size=12, bold=True)
HEAD = PatternFill("solid", fgColor="DDE6F0")
wb = Workbook()
wb.properties.creator = ""; wb.properties.lastModifiedBy = ""; wb.properties.title = "QuadSig Supplementary Data"


def cell_value(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    if isinstance(v, str):
        if v == "":
            return None
        if v in ("True", "False"):
            return v == "True"
        try:
            return int(v) if re.fullmatch(r"-?\d+", v) else float(v) if re.fullmatch(r"-?\d*\.?\d+(e-?\d+)?", v) else v
        except ValueError:
            return v
    return v


readme = [("Supplementary Data", "QuadSig: A Four-Signal, Gene-Specific Framework for Pathogenicity Assessment of p53 Missense Variants"),
          ("Contents", "One sheet per pipeline output. All values are copied unchanged from the files produced by the "
                       "pipeline (results_v2/ and data/intermediate/). Row counts are given below."),
          ("Protein", "Human p53, UniProt P04637 (isoform 1, 393 residues)"),
          ("ClinVar retrieval", clinvar_date),
          ("", "")]
for i, (name, df) in enumerate(sheets[1:], start=2):
    readme.append((f"Sheet {i}: {name}", f"{DESC[name]}  [{len(df):,} rows; columns: {', '.join(df.columns)}]"))
readme += [("", ""), ("QuadSig_Results columns", "")] + [(k, v) for k, v in RESULT_COLS.items() if k in res.columns]

ws = wb.active; ws.title = "README"
ws.column_dimensions["A"].width = 34; ws.column_dimensions["B"].width = 120
for r, (a, b) in enumerate(readme, start=1):
    ca, cb = ws.cell(r, 1, a), ws.cell(r, 2, b)
    ca.font, cb.font = (FT if r == 1 else FB), (FT if r == 1 else F)
    ca.alignment = Alignment(vertical="top"); cb.alignment = Alignment(wrap_text=True, vertical="top")

for name, df in sheets[1:]:
    ws = wb.create_sheet(name)
    for j, col in enumerate(df.columns, start=1):
        c = ws.cell(1, j, col); c.font = FB; c.fill = HEAD
        sample = [str(x) for x in df[col].head(300)]
        w = max([len(col)] + [len(s) for s in sample]) + 2
        ws.column_dimensions[get_column_letter(j)].width = min(max(w, 9), 60)
    for i, row in enumerate(df.itertuples(index=False), start=2):
        for j, v in enumerate(row, start=1):
            ws.cell(i, j, cell_value(v)).font = F
    ws.freeze_panes = "A2"
    if len(df):
        ws.auto_filter.ref = ws.dimensions
wb.save(OUT)

# ---------------- report -----------------------------------------------------------------------
print(f"\nwrote {OUT}")
for i, (name, df) in enumerate(sheets, start=1):
    print(f"  sheet {i:2d}  {name:24s} {('' if df is None else f'{len(df):,} rows'):>12s}")
cat = res.category.value_counts()
print(f"check: {len(res):,} variants; classes LP {cat.get('likely_pathogenic', 0):,} / "
      f"ambiguous {cat.get('ambiguous', 0):,} / LB {cat.get('likely_benign', 0):,}; ClinVar labelled "
      f"{len(lab)} ({(lab.clinvar_label == 'pathogenic').sum()} P, {(lab.clinvar_label == 'benign').sum()} B)")
print("problems: " + ("none" if not problems else ""))
for p in problems:
    print("  - " + p)
