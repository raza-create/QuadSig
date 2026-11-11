# QuadSig

QuadSig is a gene-specific, interpretable pipeline that scores all 7,467 possible single
amino-acid substitutions in human p53 (UniProt P04637, 393 residues). Each variant gets a
score in [0, 1], a class (likely pathogenic, ambiguous or likely benign) and the four signals
behind the score:

| Signal | Source |
|---|---|
| Population rarity | gnomAD v4 allele frequency |
| Evolutionary conservation | Ensembl orthologues, one sequence per vertebrate species (179 species) |
| Structural change | ESMFold v1 models of wild type and every mutant (C&alpha; RMSD, pLDDT) |
| Physicochemical change | hydrophobicity, volume, charge and BLOSUM62 |

A functional-domain term (DNA-binding domain, tetramerisation domain and others) is added,
and zinc-coordinating residues and gain-of-function hotspots are flagged. Classification
thresholds are chosen on ClinVar germline labels with stratified five-fold cross-validation.

This repository accompanies the manuscript *QuadSig: A Four-Signal, Gene-Specific Framework
for Pathogenicity Assessment of p53 Missense Variants* (submitted).

## Repository contents

```
code/                       pipeline steps 01-13 (sequence, mutants, ESMFold, RMSD, gnomAD,
                            orthologue alignment, conservation, physicochemical, domains,
                            zinc and gain-of-function flags) and config.py; steps 14-19 are
                            the earlier (v1) scoring and are NOT used for the paper results
recompute_conservation.py   conservation from one orthologue per vertebrate species (179)
rerun_v2.py                 QuadSig scoring, ClinVar labels, cross-validated thresholds,
                            consistency checks
make_fig23.py               Figures 2 and 3
make_supplementary.py       Additional file 1 (13-sheet Excel workbook)
data/raw/                   input sequence, orthologue alignment, annotations
data/fasta/                 mutant sequences used for structure prediction
data/intermediate/          per-variant signal tables produced by steps 01-13
results_v2/                 results reported in the paper (results.csv, cv_folds.csv,
                            parameters.json, summary.txt, ClinVar records, figures,
                            supplementary workbook)
```

The predicted structures (wild type and 7,467 mutants, about 1.5 GB) are too large for the
repository; they are attached to the GitHub release v1.0.0 as `QuadSig_structures.tar.gz`
(download it and run `tar -xzf QuadSig_structures.tar.gz` in the repository folder).

## Requirements

- Linux (tested on Ubuntu 24.04), Python 3.12
- Packages and the exact versions used: `requirements.txt`
- Structure prediction only (steps 03-04): PyTorch with CUDA 12.4, HuggingFace Transformers
  and a GPU with at least 12 GB of memory. Predicting all 7,467 mutants took 33.6 hours on
  one 12 GB consumer GPU. All other steps run without a GPU.

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Reproduce the numbers in the paper

The signal tables in `data/intermediate/` and the ClinVar records used in the paper
(`results_v2/clinvar_records.csv`) are included, so the scoring and evaluation can be
repeated without a GPU or network access:

```bash
python recompute_conservation.py      # conservation from 179 vertebrate orthologues
python rerun_v2.py --offline          # scores, classes, ROC-AUC, cross-validation, checks
python make_fig23.py
python make_supplementary.py
```

`--offline` reuses the stored ClinVar records. Without it, `rerun_v2.py` downloads the
current ClinVar records, so the evaluation numbers can differ as ClinVar changes.

## Run the full pipeline from the beginning

```bash
python code/step_01_fetch_sequence.py
python code/step_02_generate_mutants.py
python code/step_03_esmfold_wildtype.py
python code/step_04_esmfold_mutants.py      # GPU; checkpoints, so it can be resumed
python code/step_05_compute_rmsd.py
python code/step_06_genomic_coords.py
python code/step_07_gnomad_query.py
python code/step_08_download_msa.py
python code/step_09_conservation_scores.py
python code/step_10_physicochemical.py
python code/step_11_domain_annotation.py
python code/step_12_zinc_flagging.py
python code/step_13_gof_detection.py
python recompute_conservation.py            # replaces the step 09 conservation table
python rerun_v2.py
```

Steps 06-08 query public web services (Ensembl REST, gnomAD, EBI Clustal Omega), so their
output depends on the database release available at run time.

## Output

`results_v2/results.csv` has one row for each of the 7,467 variants, with the QuadSig score,
class, mechanism label, every signal value, domain, zinc and gain-of-function flags and the
ClinVar label where available.

## Data sources

UniProt (P04637), Ensembl Compara (ENSG00000141510), gnomAD v4, ClinVar and the ESMFold v1
model (huggingface.co/facebook/esmfold_v1). Please follow the terms of use of each resource.

## Archived release

Zenodo DOI: to be added after the first release.

## Licence

MIT (see `LICENSE`).
