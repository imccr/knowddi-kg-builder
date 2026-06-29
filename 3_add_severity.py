"""
Program 3 — Add DDInter severity labels to drugbank_ddi.csv
============================================================
Usage:
    python3 3_add_severity.py

Input:
    drugbank_ddi.csv            (from parse_drugbank.py)
    ddinter_downloads_code_*.csv (downloaded from ddinter.scbdd.com)

Output:
    drugbank_ddi_severity.csv   — same as drugbank_ddi.csv + severity + ddinter_id_a + ddinter_id_b
                                  severity: Minor / Moderate / Major / Unknown

Matching strategy:
    1. Normalise drug names to lowercase + apply UK/INN→US alias map
    2. Build name → DDInterID lookup from all DDInter files
    3. Join on frozenset(DDInterID_A, DDInterID_B) — order-independent
    4. Fall back to name-pair match for drugs without a DDInterID
"""

import csv, glob
from pathlib import Path
from collections import defaultdict

ROOT         = Path(__file__).parent
DDI_CSV      = ROOT / "drugbank_ddi.csv"
DDINTER_GLOB = str(ROOT / "ddinter_downloads_code_*.csv")
OUT_CSV      = ROOT / "drugbank_ddi_severity.csv"

ALIASES: dict[str, str] = {
    "aciclovir":                     "acyclovir",
    "paracetamol":                   "acetaminophen",
    "ciclosporin":                   "cyclosporine",
    "adrenaline":                    "epinephrine",
    "noradrenaline":                 "norepinephrine",
    "frusemide":                     "furosemide",
    "lignocaine":                    "lidocaine",
    "beclometasone":                 "beclomethasone",
    "colecalciferol":                "cholecalciferol",
    "glyceryl trinitrate":           "nitroglycerin",
    "suxamethonium":                 "succinylcholine",
    "pethidine":                     "meperidine",
    "hyoscine butylbromide":         "scopolamine",
    "hyoscine hydrobromide":         "scopolamine",
    "ipratropium bromide":           "ipratropium",
    "glycopyrrolate":                "glycopyrronium",
    "colistin":                      "colistimethate",
    "tenofovir disoproxil fumarate": "tenofovir",
    "p-aminosalicylic acid":         "aminosalicylic acid",
    "calcium folinate":              "leucovorin",
    "retinol":                       "vitamin a",
    "thiamine":                      "vitamin b1",
    "riboflavin":                    "vitamin b2",
    "phytomenadione":                "phytonadione",
    "l-asparaginase":                "asparaginase",
}

RANK = {"minor": 0, "moderate": 1, "major": 2}


def norm(name: str) -> str:
    return name.strip().lower()


def apply_alias(name: str) -> str:
    n = norm(name)
    return ALIASES.get(n, n)


def load_ddinter(pattern: str):
    """
    Returns:
        name_to_id  : normed_drug_name → DDInterID  (e.g. "abacavir" → "DDInter1")
        id_pair_map : frozenset({DDInterID_A, DDInterID_B}) → severity level
        name_pair_map: frozenset({normed_name_a, normed_name_b}) → severity level (fallback)
    """
    name_to_id:   dict[str, str]         = {}
    id_pair_map:  dict[frozenset, str]   = {}
    name_pair_map: dict[frozenset, str]  = {}

    for fpath in sorted(glob.glob(pattern)):
        with open(fpath, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                id_a   = row["DDInterID_A"].strip()
                name_a = apply_alias(row["Drug_A"])
                id_b   = row["DDInterID_B"].strip()
                name_b = apply_alias(row["Drug_B"])
                level  = row["Level"].strip().capitalize()

                name_to_id.setdefault(name_a, id_a)
                name_to_id.setdefault(name_b, id_b)

                id_key   = frozenset({id_a, id_b})
                name_key = frozenset({name_a, name_b})

                for mapping, key in ((id_pair_map, id_key), (name_pair_map, name_key)):
                    existing = mapping.get(key)
                    if existing is None or RANK.get(level.lower(), -1) > RANK.get(existing.lower(), -1):
                        mapping[key] = level

    return name_to_id, id_pair_map, name_pair_map


def main():
    print("Loading DDInter files …")
    name_to_id, id_pair_map, name_pair_map = load_ddinter(DDINTER_GLOB)
    print(f"  {len(name_to_id):,} unique drug names → DDInterID mappings")
    print(f"  {len(id_pair_map):,} unique DDInterID pairs with severity")

    matched_by_id = matched_by_name = unmatched = 0
    rows_out = []

    with open(DDI_CSV, encoding="utf-8") as f:
        reader    = csv.DictReader(f)
        fieldnames = reader.fieldnames + ["ddinter_id_a", "ddinter_id_b", "severity"]

        for row in reader:
            name_a = apply_alias(row["drug_name"])
            name_b = apply_alias(row["partner_name"])

            id_a = name_to_id.get(name_a, "")
            id_b = name_to_id.get(name_b, "")

            level = "Unknown"

            # Primary: match by DDInter ID pair
            if id_a and id_b:
                level = id_pair_map.get(frozenset({id_a, id_b}), "Unknown")
                if level != "Unknown":
                    matched_by_id += 1

            # Fallback: match by normalised name pair
            if level == "Unknown":
                level = name_pair_map.get(frozenset({name_a, name_b}), "Unknown")
                if level != "Unknown":
                    matched_by_name += 1

            if level == "Unknown":
                unmatched += 1

            row["ddinter_id_a"] = id_a
            row["ddinter_id_b"] = id_b
            row["severity"]     = level
            rows_out.append(row)

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    total = len(rows_out)
    print(f"\nResults:")
    print(f"  Total DDI rows        : {total:,}")
    print(f"  Matched by DDInterID  : {matched_by_id:,} ({100*matched_by_id//total}%)")
    print(f"  Matched by name       : {matched_by_name:,} ({100*matched_by_name//total}%)")
    print(f"  Unknown               : {unmatched:,} ({100*unmatched//total}%)")
    print(f"\nOutput → {OUT_CSV}")

    from collections import Counter
    dist = Counter(r["severity"] for r in rows_out)
    print("\nSeverity distribution:")
    for level, count in sorted(dist.items(), key=lambda x: -x[1]):
        print(f"  {level:<10} {count:>8,}")


if __name__ == "__main__":
    main()
