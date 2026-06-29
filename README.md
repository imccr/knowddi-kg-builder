# KnowDDI-Nepal

Nepal-specific extension of [KnowDDI](https://github.com/xzenglab/KnowDDI) — a Graph Neural Network framework for Drug-Drug Interaction (DDI) prediction.

This repository contains:
- The DrugBank knowledge graph extended with Nepal NLEM drugs
- Scripts to add new Nepal drugs to the DDI and BKG graph at any time
- The original KnowDDI training pipeline

The trained model and web application live in the companion repo: **[MediSafe-Nepal](https://github.com/your-username/medisafe-nepal)**.

---

## What is KnowDDI?

KnowDDI combines a Drug-Drug Interaction graph with a Biological Knowledge Graph (BKG) from [Hetionet](https://het.io/) to train a GNN that predicts 86 types of drug interactions. This repo extends the original DrugBank dataset of **1,710 drugs** to cover drugs on Nepal's National List of Essential Medicines (NLEM 2021).

---

## Nepal Drug Coverage

| Category | Count | Details |
|---|---|---|
| Drugs already in KnowDDI | 286 | Matched from NLEM via name/synonym — see `MAPPED_286.json` |
| New drugs appended | 34 | Node IDs 34124–34157 — see `APPENDED_34.json` |
| Unmapped (no DrugBank match) | 31 | See `UNMAPPED_31.json` for details and reasons |
| Total NLEM coverage | 320 / 324 | |

---

## Repository Structure

```
knowddi-nepal/
│
├── 1_find_drug_ids.py        # Step 1: find DrugBank IDs for new drug names
├── 2_pull_edges.py           # Step 2: extract DDI + BKG edges for new drugs
│
├── MAPPED_286.json           # 286 NLEM drugs already in original KnowDDI
├── APPENDED_34.json           # 34 new drugs added (format: name → {drugbank_id, knowddi_id, ...})
├── UNMAPPED_31.json          # 31 drugs initially unmapped (now resolved in APPENDED_7)
├── nepal_drugs_extracted.json # All 324 NLEM 2021 drug names
│
├── drugbank_drugs.csv        # Parsed from DrugBank XML: id, name, type, synonyms
├── drugbank_targets.csv      # Parsed from DrugBank XML: drug → gene targets
│   # drugbank_ddi.csv        # 87 MB — gitignored, regenerate with parse_drugbank.py
│   # drugbank.xml.zip        # 85 MB — gitignored, download from drugbank.ca
│
├── data/
│   ├── drugbank/             # DrugBank graph files used for training
│   │   ├── train.txt         # DDI edges: head_id  tail_id  rel_type_id
│   │   ├── valid.txt
│   │   ├── test.txt
│   │   └── BKG_file.txt      # Biological KG edges: entity_id  entity_id  rel_id
│   └── BioSNAP/              # BioSNAP benchmark dataset (unchanged)
│
├── raw_data/
│   ├── Drugbank/             # Node/entity/relation mappings
│   │   ├── node2id.json      # DrugBank name → integer node ID (0–34157)
│   │   ├── BKG_entity2Id.json # Hetionet entity → integer ID (1710–34123)
│   │   ├── id2rel.txt        # 86 DDI interaction templates
│   │   └── relation_type_drug.json  # 23 BKG relation type names
│   ├── BioSNAP/
│   └── hetionet/
│       └── hetionet-v1.0-nodes.tsv  # Gene symbol → Entrez Gene ID
│
├── pytorch/                  # Original KnowDDI training code (upstream)
│
├── MediSafe_Dataset_Expansion_Guide.pdf
├── requirements.txt
└── .gitignore
```

---

## Workflow: Adding New Nepal Drugs

Use this workflow any time you discover new drugs sold in Nepal that need to be added to the interaction graph.

### Prerequisites

```bash
pip install -r requirements.txt
# Also required: drugbank_drugs.csv and drugbank_targets.csv
# If missing, regenerate them:
#   unzip drugbank.xml.zip
#   python3 parse_drugbank.py
```

### Step 1 — Find DrugBank IDs

Create a plain text file with your new drug names (one per line):

```
# new_drugs.txt
clonazepam
verapamil
amiodarone
```

Run the finder, skipping drugs already mapped:

```bash
python3 1_find_drug_ids.py \
    --file new_drugs.txt \
    --skip MAPPED_286.json APPENDED_34.json

# For fuzzy/uncertain matches, confirm interactively:
python3 1_find_drug_ids.py \
    --file new_drugs.txt \
    --skip MAPPED_286.json APPENDED_34.json \
    --interactive
```

**You can also pass `nepal_drugs_extracted.json` directly** — the `--skip` flag automatically excludes already-mapped drugs:

```bash
python3 1_find_drug_ids.py \
    --file nepal_drugs_extracted.json \
    --skip MAPPED_286.json APPENDED_34.json
```

Output → `new_drug_ids.json` (same schema as `APPENDED_34.json`)

### Step 2 — Extract DDI + BKG Edges

```bash
python3 2_pull_edges.py --input new_drug_ids.json
```

Output — three files ready to copy-paste:

| File | Paste into |
|---|---|
| `out_node2id_entries.json` | `raw_data/Drugbank/node2id.json` |
| `out_train_edges.txt` | `data/drugbank/train.txt` |
| `out_BKG_edges.txt` | `data/drugbank/BKG_file.txt` |

### Step 3 — Update Mapping Files

Merge `new_drug_ids.json` into `APPENDED_34.json` so the `--skip` list stays current for the next batch.

### Step 4 — Retrain

See **[MediSafe-Nepal](https://github.com/your-username/medisafe-nepal)** for the training pipeline. The updated `data/drugbank/` files are the direct input to `train.py`.

---

## Match Tiers (how Program 1 finds drugs)

| Tier | Method | Example |
|---|---|---|
| `exact_name` | Exact match on DrugBank canonical name | `acetylcysteine` → Acetylcysteine |
| `exact_synonym` | Exact match on any synonym | `dihydroartemisinin` → Artenimol |
| `substr_name` | Query is prefix/substring of name | `atracurium` → Atracurium besylate |
| `substr_synonym` | Query substring of a synonym | `beclometasone` → Beclomethasone dipropionate |
| `fuzzy` ⚠ | Character similarity ≥ 0.78 | `fomipezole` → Fomepizole |
| NOT FOUND | Absent from DrugBank 2017 | `mesna`, `piperaquine` |

Fuzzy matches are flagged `NEEDS_REVIEW: true` in the JSON output. Run with `--interactive` to confirm them manually.

---

## Training

The `pytorch/` folder contains the original upstream KnowDDI training code. The refactored version used in production is in **MediSafe-Nepal**.

To train with the original code:

```bash
cd pytorch
python3 train.py
```

---

## Data Sources

| Source | Version | License |
|---|---|---|
| [DrugBank](https://www.drugbank.ca) | 5.0 (2017-12-20) | Academic/Non-commercial |
| [Hetionet](https://het.io/) | v1.0 | CC0 |
| [BioSNAP](http://snap.stanford.edu/biodata/) | — | Academic |
| [DDInter](https://ddinter.scbdd.com) | 2.0 | Academic |
| Nepal NLEM | 2021 | Government of Nepal |

> **DrugBank XML is not included** in this repo due to licensing. Download `drugbank.xml` from [drugbank.ca](https://www.drugbank.ca/releases/latest) (free academic registration required) and place it in the root directory before running `parse_drugbank.py`.

> **DDInter CSVs are not included** in this repo. Download all category files from [ddinter.scbdd.com/download](https://ddinter.scbdd.com/download/) and place them inside a `ddinter_data/` folder in the root directory before running `3_add_severity.py`.

---

## Citation



```bibtex
@article{knowddi2022,
  title   = {KnowDDI: ...},
  author  = {...},
  journal = {...},
  year    = {2022}
}
```

---

## Related

- [MediSafe-Nepal](https://github.com/your-username/medisafe-nepal) — FastAPI + Streamlit web application for DDI prediction
