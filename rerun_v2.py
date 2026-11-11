#!/usr/bin/env python3
"""QuadSig v2 rerun: scoring + ClinVar consensus labels + cross-validated thresholds +
seven consistency checks + Figures 2 and 3, all from the existing signal files.
Design (2026-11-09): base formula of step 14 with NO minimum scores and NO extra rules;
ClinVar consensus rule; Youden cut-off with stratified 5-fold CV; ambiguous band = range of
the five fold cut-offs. Reads data/intermediate/* (unchanged). Writes ONLY to results_v2/.
Usage: python3 rerun_v2.py   |   python3 rerun_v2.py --offline (reuse saved ClinVar records)
"""
import argparse, csv, datetime, json, math, os, sys, time, urllib.parse, urllib.request
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = os.path.dirname(os.path.abspath(__file__))
INTER = os.path.join(ROOT, "data", "intermediate")
OUT = os.path.join(ROOT, "results_v2")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

W = dict(pop=0.15, cons=0.25, rmsd=0.20, plddt=0.10, phys=0.15, inter=0.075)
DOMAIN_ADD = {"DBD": 0.15, "TET": 0.08, "TAD": 0.02, "PRR": 0.0, "CTD": -0.05, "LINKER": -0.03}
HOTSPOTS = ["R175H", "G245S", "R248W", "R248Q", "R249S", "R273H", "R273C"]
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
B62_AA = "ARNDCQEGHILKMFPSTWYV"
B62 = {
    "A": [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    "R": [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    "N": [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    "D": [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    "C": [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    "Q": [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    "E": [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    "G": [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    "H": [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    "I": [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    "L": [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    "K": [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    "M": [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    "F": [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    "P": [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    "S": [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    "T": [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    "W": [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    "Y": [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    "V": [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
}


def blosum(wt, mt):
    return B62[wt][B62_AA.index(mt)]


def load_signals(say):
    rd = lambda f: pd.read_csv(os.path.join(INTER, f), dtype=str)
    muts = [l.strip() for l in open(os.path.join(INTER, "mutations.txt")) if l.strip()]
    df = pd.DataFrame({"mutation": muts})
    df["position"] = df.mutation.str[1:-1].astype(int)
    df["wt_aa"], df["mut_aa"] = df.mutation.str[0], df.mutation.str[-1]
    st = rd("structural_scores.csv")[["mutation", "rmsd", "plddt_wt", "plddt_mut", "plddt_change"]]
    po = rd("population_scores.csv")[["mutation", "population_score", "allele_frequency"]]
    co = rd("conservation_scores.csv")[["position", "conservation_score"]]
    co["position"] = co.position.astype(int)
    ph = rd("physicochemical_scores.csv")[["mutation", "physicochemical_score"]]
    do = rd("domain_annotations.csv")[["mutation", "domain_code", "domain_name"]]
    zn = rd("zinc_flags.csv")[["mutation", "is_zinc_coord"]]
    gf = rd("gof_flags.csv")
    gf = gf[[c for c in ["mutation", "is_gof", "gof_mechanism", "structural_divergence_note"] if c in gf]]
    for t in [st, po, ph, do, zn, gf]:
        df = df.merge(t, on="mutation", how="left")
    df = df.merge(co, on="position", how="left")
    defaults = {"population_score": 1.0, "conservation_score": 0.5, "rmsd": 0.0,
                "plddt_change": 0.0, "physicochemical_score": 0.0, "allele_frequency": 0.0}
    for c, v in defaults.items():
        n_missing = int(pd.to_numeric(df[c], errors="coerce").isna().sum())
        if n_missing:
            say(f"  note: {n_missing} missing values in {c} -> default {v} (as in original step 14)")
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(v)
    df["domain_code"] = df.domain_code.fillna("LINKER")
    for c in ["is_zinc_coord", "is_gof"]:
        df[c] = df[c].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    return df


def quadsig(df):
    rmsd_norm = 1.0 / (1.0 + np.exp(-1.5 * (df.rmsd - 2.0)))
    pen = np.where(df.plddt_change < 0, np.minimum(df.plddt_change.abs() / 30.0, 1.0), 0.0)
    s = (W["pop"] * df.population_score + W["cons"] * df.conservation_score + W["rmsd"] * rmsd_norm
         + W["plddt"] * pen + W["phys"] * df.physicochemical_score
         + W["inter"] * df.conservation_score * df.physicochemical_score)
    s = s + df.domain_code.map(DOMAIN_ADD).fillna(0.0)
    af = df.allele_frequency
    s = np.where(af > 0.01, s * 0.2, np.where(af > 0.001, s * 0.4, np.where(af > 0.0001, s * 0.6, s)))
    return np.round(np.clip(s, 0, 1), 4), rmsd_norm, pen


def mechanism(r):
    p = []
    if r.is_zinc_coord: p.append(f"Zinc-coordinating residue {r.wt_aa}{r.position}")
    if r.is_gof: p.append("Curated GOF hotspot")
    if r.domain_code == "DBD": p.append("DNA-binding domain")
    if r.conservation_score > 0.8: p.append(f"Highly conserved position (cons={r.conservation_score:.2f})")
    if r.rmsd > 2.0: p.append(f"Large predicted structural change (RMSD={r.rmsd:.2f})")
    if r.plddt_change < -10: p.append(f"Confidence drop (pLDDT change={r.plddt_change:.1f})")
    if r.physicochemical_score > 0.5: p.append("Chemically drastic substitution")
    if r.population_score < 0.5: p.append("Observed in gnomAD at AF >= 1e-3")
    return "; ".join(p)


def get_json(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "QuadSig/2.0"}), timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 + 2 * k)


def record_label(desc):
    d = desc.lower()
    if "conflicting" in d:
        return "conflicting"
    if "pathogenic" in d and "benign" not in d and "uncertain" not in d:
        return "P"
    if "benign" in d and "pathogenic" not in d and "uncertain" not in d:
        return "B"
    return "U"


def parse_pchange(title):
    if "(p." not in title:
        return ""
    pc = title.split("(p.")[1].split(")")[0]
    wt, mt, pos = AA3.get(pc[:3], ""), AA3.get(pc[-3:], ""), pc[3:-3]
    return wt + pos + mt if (wt and mt and pos.isdigit()) else ""


def fetch_clinvar(say):
    term = "TP53[gene] AND single_nucleotide_variant[Type of variation]"
    url = f"{EUTILS}/esearch.fcgi?db=clinvar&term={urllib.parse.quote(term)}&retmax=20000&retmode=json"
    ids = get_json(url)["esearchresult"]["idlist"]
    say(f"  ClinVar TP53 SNV records found: {len(ids)}")
    rows = []
    for i in range(0, len(ids), 100):
        res = get_json(f"{EUTILS}/esummary.fcgi?db=clinvar&id={','.join(ids[i:i+100])}&retmode=json")["result"]
        for vid in ids[i:i + 100]:
            e = res.get(vid, {})
            genes = [g.get("symbol", "") for g in e.get("genes", [])]
            g = e.get("germline_classification", {}) or {}
            rows.append(dict(uid=vid, title=e.get("title", ""), genes=";".join(genes),
                             germline=g.get("description", ""), review=g.get("review_status", ""),
                             last_evaluated=g.get("last_evaluated", "")))
        time.sleep(0.4)
    return pd.DataFrame(rows)


def consensus(records):
    r = records[(records.genes.str.split(";").str[0] == "TP53") & (records.germline != "")].copy()
    r["mutation"] = r.title.map(parse_pchange)
    r = r[r.mutation != ""]
    r["rec"] = r.germline.map(record_label)
    out = []
    for m, g in r.groupby("mutation"):
        s = set(g.rec)
        lab = ("pathogenic" if ("P" in s and "B" not in s) else
               "benign" if ("B" in s and "P" not in s) else "uncertain")
        out.append(dict(mutation=m, clinvar_label=lab, n_records=len(g),
                        clinvar_significance=" | ".join(sorted(set(g.germline)))))
    return pd.DataFrame(out), len(r)


def youden(scores, y):
    best_t, best_j = 0.5, -1
    for t in np.unique(scores):
        pred = scores >= t
        sens = (pred & (y == 1)).sum() / max((y == 1).sum(), 1)
        spec = (~pred & (y == 0)).sum() / max((y == 0).sum(), 1)
        if sens + spec - 1 > best_j:
            best_j, best_t = sens + spec - 1, t
    return float(best_t)


def metrics(pred, y):
    tp, fn = int((pred & (y == 1)).sum()), int((~pred & (y == 1)).sum())
    fp, tn = int((pred & (y == 0)).sum()), int((~pred & (y == 0)).sum())
    return dict(TP=tp, FN=fn, FP=fp, TN=tn, accuracy=(tp + tn) / len(y), sensitivity=tp / max(tp + fn, 1),
                specificity=tn / max(tn + fp, 1),
                balanced_accuracy=0.5 * (tp / max(tp + fn, 1) + tn / max(tn + fp, 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="reuse results_v2/clinvar_records.csv")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    log = open(os.path.join(OUT, "summary.txt"), "w")

    def say(*s):
        m = " ".join(str(x) for x in s); print(m, flush=True); log.write(m + "\n"); log.flush()

    say(f"QuadSig v2 rerun  {datetime.datetime.now():%Y-%m-%d %H:%M}")
    say("\n[1] Scoring (no minimum scores, no extra rules)")
    df = load_signals(say)
    df["quadsig_score"], df["rmsd_norm"], df["plddt_penalty"] = quadsig(df)
    say(f"  {len(df)} mutations; score range {df.quadsig_score.min():.3f}-{df.quadsig_score.max():.3f}, "
        f"mean {df.quadsig_score.mean():.3f}")

    say("\n[2] ClinVar (germline classification, consensus rule)")
    rec_path = os.path.join(OUT, "clinvar_records.csv")
    if a.offline and os.path.exists(rec_path):
        records = pd.read_csv(rec_path, dtype=str).fillna("")
        say(f"  reused {rec_path}")
    else:
        records = fetch_clinvar(say)
        records.to_csv(rec_path, index=False)
        say(f"  fetched {datetime.date.today()}; saved {rec_path}")
    cv, n_missense_records = consensus(records)
    df = df.merge(cv, on="mutation", how="left")
    df["clinvar_label"] = df.clinvar_label.fillna("")
    lab = df[df.clinvar_label.isin(["pathogenic", "benign"])].copy()
    y = (lab.clinvar_label == "pathogenic").astype(int).to_numpy()
    s = lab.quadsig_score.to_numpy()
    say(f"  missense records with germline classification: {n_missense_records}; "
        f"unique protein changes: {len(cv)} (pathogenic {int((cv.clinvar_label=='pathogenic').sum())}, "
        f"benign {int((cv.clinvar_label=='benign').sum())}, uncertain {int((cv.clinvar_label=='uncertain').sum())})")
    say(f"  labelled variants matched to the 7,467 set: {len(lab)} ({int(y.sum())} pathogenic, {int((1-y).sum())} benign)")

    say("\n[3] Threshold and performance (Youden cut-off; stratified 5-fold CV, seed 42)")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof = np.zeros(len(y), dtype=bool); fold_t, fold_rows = [], []
    for k, (tr, te) in enumerate(skf.split(s, y), 1):
        t = youden(s[tr], y[tr]); fold_t.append(t)
        oof[te] = s[te] >= t
        m = metrics(s[te] >= t, y[te]); fold_rows.append(dict(fold=k, threshold=t, **m))
    cvm = metrics(oof, y)
    t_all = youden(s, y); ins = metrics(s >= t_all, y)
    auc = roc_auc_score(y, s)
    pd.DataFrame(fold_rows).to_csv(os.path.join(OUT, "cv_folds.csv"), index=False)
    lo, hi = min(fold_t), max(fold_t)
    say(f"  ROC-AUC (threshold-free): {auc:.3f}")
    say(f"  fold thresholds: {', '.join(f'{t:.4f}' for t in fold_t)}")
    say(f"  CROSS-VALIDATED (out-of-fold, n={len(y)}): accuracy {100*cvm['accuracy']:.1f}%, "
        f"sensitivity {100*cvm['sensitivity']:.1f}%, specificity {100*cvm['specificity']:.1f}%, "
        f"balanced accuracy {100*cvm['balanced_accuracy']:.1f}%  [TP {cvm['TP']} FN {cvm['FN']} FP {cvm['FP']} TN {cvm['TN']}]")
    say(f"  IN-SAMPLE (threshold {t_all:.4f} fit on all labelled): accuracy {100*ins['accuracy']:.1f}%, "
        f"sensitivity {100*ins['sensitivity']:.1f}%, specificity {100*ins['specificity']:.1f}%")
    say(f"  Final classification: likely pathogenic >= {hi:.4f}; ambiguous {lo:.4f}-{hi:.4f} "
        f"(range of fold thresholds); likely benign < {lo:.4f}")
    df["category"] = np.where(df.quadsig_score >= hi, "likely_pathogenic",
                              np.where(df.quadsig_score < lo, "likely_benign", "ambiguous"))
    df["clinvar_agreement"] = ""
    lp, lb = df.category == "likely_pathogenic", df.category == "likely_benign"
    df.loc[(df.clinvar_label == "pathogenic") & lp | (df.clinvar_label == "benign") & lb, "clinvar_agreement"] = "AGREE"
    df.loc[df.clinvar_label.isin(["pathogenic", "benign"]) & (df.clinvar_agreement == ""), "clinvar_agreement"] = "DISAGREE"
    cats = df.category.value_counts()
    say("  categories: " + ", ".join(f"{c} {cats.get(c, 0)} ({100*cats.get(c, 0)/len(df):.1f}%)"
                                     for c in ["likely_pathogenic", "ambiguous", "likely_benign"]))
    vus = df[df.clinvar_label == "uncertain"]
    say(f"  ClinVar uncertain/conflicting protein changes in the set: {len(vus)} -> "
        + ", ".join(f"{c} {int((vus.category == c).sum())}" for c in ["likely_pathogenic", "ambiguous", "likely_benign"]))

    say("\n[4] Seven consistency checks (targets set before the rerun; reported as met / not met)")
    df["blosum62"] = [blosum(w, m) for w, m in zip(df.wt_aa, df.mut_aa)]
    rho_c = spearmanr(df.conservation_score, df.quadsig_score).correlation
    rho_b = spearmanr(df.blosum62, df.quadsig_score).correlation
    dom_mean = df.groupby("domain_code").quadsig_score.mean().sort_values(ascending=False)
    zn_mean = df[df.is_zinc_coord].quadsig_score.mean()
    hs = df.set_index("mutation").loc[HOTSPOTS]
    n_hs = int((hs.category == "likely_pathogenic").sum())
    dbd = df.domain_code == "DBD"
    p_dbd = 100 * (df[dbd].category == "likely_pathogenic").mean()
    p_non = 100 * (df[~dbd].category == "likely_pathogenic").mean()
    checks = [("ClinVar agreement, cross-validated", f"{100*cvm['accuracy']:.1f}%", ">= 85%", cvm["accuracy"] >= 0.85),
              ("Conservation correlation (Spearman)", f"{rho_c:.2f}", ">= 0.50", rho_c >= 0.5),
              ("Domain gradient (DBD highest)", f"{dom_mean.index[0]} = {dom_mean.iloc[0]:.2f}", "DBD highest", dom_mean.index[0] == "DBD"),
              ("Zinc-residue mean score", f"{zn_mean:.2f}", ">= 0.85", zn_mean >= 0.85),
              ("GOF hotspots likely pathogenic", f"{n_hs}/7", "7/7", n_hs == 7),
              ("BLOSUM62 anti-correlation (Spearman)", f"{rho_b:.2f}", "<= -0.30", rho_b <= -0.30),
              ("DBD pathogenic enrichment", f"{p_dbd:.0f}% vs {p_non:.0f}%", "DBD > non-DBD", p_dbd > p_non)]
    for c, r, t, ok in checks:
        say(f"  {c:40s} {r:>16s}  target {t:12s} {'MET' if ok else 'NOT MET'}")
    say("  domain means: " + ", ".join(f"{d} {v:.2f} (n={int((df.domain_code == d).sum())})" for d, v in dom_mean.items()))
    say(f"  zinc substitutions: n={int(df.is_zinc_coord.sum())}, min {df[df.is_zinc_coord].quadsig_score.min():.3f}, "
        f"{int((df[df.is_zinc_coord].category == 'likely_pathogenic').sum())} likely pathogenic")
    say("  hotspots: " + ", ".join(f"{m} {hs.loc[m, 'quadsig_score']:.2f} ({hs.loc[m, 'category'].replace('likely_', '')}; "
                                   f"ClinVar {hs.loc[m, 'clinvar_label'] or 'none'}; RMSD {hs.loc[m, 'rmsd']:.2f})" for m in HOTSPOTS))
    p72 = df.set_index("mutation").loc["P72R"]
    say(f"  P72R: score {p72.quadsig_score:.3f}, {p72.category}, ClinVar {p72.clinvar_label or 'none'}, AF {p72.allele_frequency}")
    say(f"  gnomAD: {int((df.allele_frequency > 0).sum())} of {len(df)} observed; hotspots observed: "
        + (", ".join(m for m in HOTSPOTS if hs.loc[m, 'allele_frequency'] > 0) or "none"))

    df["mechanism"] = df.apply(mechanism, axis=1)
    cols = ["mutation", "position", "wt_aa", "mut_aa", "quadsig_score", "category", "population_score",
            "allele_frequency", "conservation_score", "rmsd", "rmsd_norm", "plddt_wt", "plddt_mut",
            "plddt_change", "plddt_penalty", "physicochemical_score", "blosum62", "domain_code",
            "domain_name", "is_zinc_coord", "is_gof", "mechanism", "clinvar_label", "n_records",
            "clinvar_significance", "clinvar_agreement"]
    df[cols].to_csv(os.path.join(OUT, "results.csv"), index=False)
    json.dump(dict(weights=W, domain_add=DOMAIN_ADD, threshold_low=lo, threshold_high=hi,
                   youden_all=t_all, fold_thresholds=fold_t, auc=auc, cv=cvm, in_sample=ins),
              open(os.path.join(OUT, "parameters.json"), "w"), indent=1)

    say("\n[5] Figures from results_v2/results.csv")
    figures(df, lo, hi, rho_c, rho_b, cvm, say)
    say(f"\nDone. Outputs in {OUT}")


def figures(df, lo, hi, rho_c, rho_b, cvm, say):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 600, "savefig.bbox": "tight"})
    CAT = [("likely_benign", "Likely benign", "#3b6fb6"), ("ambiguous", "Ambiguous", "#9a9a9a"),
           ("likely_pathogenic", "Likely pathogenic", "#c8413b")]
    DOMC = {"DBD": "#3b6fb6", "TET": "#e07b39", "LINKER": "#4c9a52", "CTD": "#8e5ea2", "TAD": "#d4a017", "PRR": "#2a9d8f"}
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 0.55], hspace=0.55, wspace=0.28)
    ax = fig.add_subplot(gs[0, 0]); bins = np.arange(0, 1.0001, 0.02)
    ax.hist([df[df.category == c].quadsig_score for c, _, _ in CAT], bins=bins, stacked=True,
            color=[c for _, _, c in CAT], label=[f"{l} (n={int((df.category == c).sum()):,})" for c, l, _ in CAT],
            edgecolor="white", linewidth=0.3)
    for t in sorted({round(lo, 4), round(hi, 4)}):
        ax.axvline(t, ls="--", lw=0.8, color="k")
    ax.set_xlabel("QuadSig score"); ax.set_ylabel("Number of variants"); ax.set_xlim(0, 1)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
    ax = fig.add_subplot(gs[0, 1])
    order = df.groupby("domain_code").quadsig_score.mean().sort_values(ascending=False).index.tolist()
    data = [df[df.domain_code == d].quadsig_score.to_numpy() for d in order]
    vp = ax.violinplot(data, showextrema=False, widths=0.8)
    for b, d in zip(vp["bodies"], order):
        b.set_facecolor(DOMC.get(d, "#777")); b.set_alpha(0.35); b.set_edgecolor(DOMC.get(d, "#777"))
    ax.boxplot(data, widths=0.18, showfliers=False, medianprops=dict(color="k", lw=1),
               boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for i, d in enumerate(data, 1):
        ax.text(i, 1.03, f"{d.mean():.2f}", ha="center", fontsize=7, fontweight="bold")
    ax.axhline(hi, ls="--", lw=0.8, color="k")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([f"{d}\nn={len(x):,}" for d, x in zip(order, data)], fontsize=6.5)
    ax.set_ylim(0, 1.1); ax.set_ylabel("QuadSig score")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    ax = fig.add_subplot(gs[1, :]); left = 0
    for c, l, col in [CAT[2], CAT[1], CAT[0]]:
        n = int((df.category == c).sum()); w = 100 * n / len(df)
        ax.barh(0, w, left=left, color=col, edgecolor="white")
        if w > 4:
            ax.text(left + w / 2, 0, f"{l}\n{n:,} ({w:.1f}%)", ha="center", va="center", color="white", fontsize=7, fontweight="bold")
        left += w
    ax.set_xlim(0, 100); ax.set_yticks([]); ax.set_xlabel("Percentage of all 7,467 variants")
    ax.spines["left"].set_visible(False)
    ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
    for e in ["png", "pdf"]:
        fig.savefig(os.path.join(OUT, f"fig2.{e}"))
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    ax = axes[0]
    for d in order:
        x = df[df.domain_code == d]
        ax.scatter(x.conservation_score, x.quadsig_score, s=2, alpha=0.35, color=DOMC.get(d, "#777"), lw=0, label=d, rasterized=True)
    ax.set_xlabel("Conservation score"); ax.set_ylabel("QuadSig score"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=5.5, markerscale=3, loc="lower right", ncol=2)
    ax.set_title(f"a   Spearman $\\rho$ = {rho_c:.2f}", loc="left", fontsize=8, fontweight="bold")
    ax = axes[1]
    vals = sorted(df.blosum62.unique())
    bx = ax.boxplot([df[df.blosum62 == v].quadsig_score for v in vals], showfliers=False, patch_artist=True,
                    medianprops=dict(color="k"), widths=0.6)
    cm = plt.get_cmap("RdYlGn")
    for i, p in enumerate(bx["boxes"]):
        p.set_facecolor(cm(i / max(len(vals) - 1, 1))); p.set_linewidth(0.6)
    ax.set_xticks(range(1, len(vals) + 1)); ax.set_xticklabels([str(v) for v in vals])
    ax.set_xlabel("BLOSUM62 substitution score"); ax.set_ylabel("QuadSig score"); ax.set_ylim(0, 1)
    ax.set_title(f"b   Spearman $\\rho$ = {rho_b:.2f}", loc="left", fontsize=8, fontweight="bold")
    ax = axes[2]; bins = np.arange(0, 1.0001, 0.025)
    for lab, col in [("benign", "#4caf7a"), ("pathogenic", "#d9534f")]:
        x = df[df.clinvar_label == lab].quadsig_score
        ax.hist(x, bins=bins, alpha=0.6, color=col, label=f"ClinVar {lab} (n={len(x)})")
    ax.axvline(hi, ls="--", lw=0.8, color="k", label=f"Threshold {hi:.2f}")
    ax.set_xlabel("QuadSig score"); ax.set_ylabel("Number of variants"); ax.set_xlim(0, 1)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.4)
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")
    ax.set_title(f"c   CV accuracy {100*cvm['accuracy']:.1f}%", loc="left", fontsize=8, fontweight="bold")
    fig.tight_layout()
    for e in ["png", "pdf"]:
        fig.savefig(os.path.join(OUT, f"fig3.{e}"))
    plt.close(fig)
    say("  wrote fig2.png/.pdf and fig3.png/.pdf (600 dpi)")


if __name__ == "__main__":
    main()
