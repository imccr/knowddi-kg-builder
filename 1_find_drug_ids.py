"""
Program 1 — Find DrugBank IDs for new Nepal drug names
========================================================
Usage:
    python3 1_find_drug_ids.py "drug one" "drug two" ...
    python3 1_find_drug_ids.py --file drug_names.txt          # plain text, one name per line
    python3 1_find_drug_ids.py --file nepal_drugs_extracted.json   # JSON list also supported
    python3 1_find_drug_ids.py --file nepal_drugs_extracted.json --skip MAPPED_286.json APPENDED_7.json
    python3 1_find_drug_ids.py --file list.txt --interactive   # confirm fuzzy hits

Typical workflow for Nepal NLEM expansion:
    python3 1_find_drug_ids.py \\
        --file nepal_drugs_extracted.json \\
        --skip MAPPED_286.json APPENDED_7.json \\
        --interactive

No hardcoded alias tables.  Matching tiers (applied in order, stops at first hit):
  1. exact_name      — norm(query) == norm(drugbank canonical name)
  2. exact_synonym   — norm(query) == norm(any synonym)
  3. substr_name     — query tokens all contained in drugbank name (or vice versa)
  4. substr_synonym  — same but against synonyms
  5. token_overlap   — ≥ 80% of query tokens found in any name/synonym
  6. fuzzy           — difflib similarity ≥ 0.82 against any name or synonym

Tiers 1–4 are auto-accepted.
Tier 5–6 show the top candidates; with --interactive you confirm each one.
Without --interactive they are written with match="fuzzy" and marked ⚠ in console.

Output: new_drug_ids.json  (same schema as APPENDED_7.json)
"""

import csv, re, json, sys, argparse
from pathlib import Path
from difflib import SequenceMatcher

ROOT = Path(__file__).parent

DRUGS_CSV     = ROOT / "drugbank_drugs.csv"
APPENDED_JSON = ROOT / "APPENDED_7.json"
NODE2ID_JSON  = ROOT / "raw_data/Drugbank/node2id.json"
OUT_JSON      = ROOT / "new_drug_ids.json"

FUZZY_THRESHOLD   = 0.78   # minimum similarity ratio for fuzzy tier
TOKEN_MIN_OVERLAP = 0.80   # fraction of query tokens that must appear in candidate


# ── normalisation ─────────────────────────────────────────────────────────────

def norm(s: str) -> str:
    """Lowercase, collapse every non-alphanumeric run to a single space."""
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()

def tokenize(s: str) -> set[str]:
    """Return set of meaningful tokens (length ≥ 4, avoids element symbols and abbrevs)."""
    return {t for t in norm(s).split() if len(t) >= 4}


# ── index building ────────────────────────────────────────────────────────────

def load_index() -> tuple[dict, dict, list]:
    """
    Returns
      name_idx   : norm_name  → [(db_id, canonical)]
      syn_idx    : norm_syn   → [(db_id, canonical)]
      flat_list  : [(norm_text, db_id, canonical, is_name)]
                   flat list for token-overlap and fuzzy search
    """
    name_idx: dict[str, list] = {}
    syn_idx:  dict[str, list] = {}
    flat_list: list[tuple[str, str, str, bool]] = []

    with open(DRUGS_CSV, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            db_id = row["drugbank_id"]
            cname = row["name"]
            nn    = norm(cname)
            name_idx.setdefault(nn, []).append((db_id, cname))
            flat_list.append((nn, db_id, cname, True))

            for syn in row["synonyms"].split("|"):
                ns = norm(syn)
                if ns and ns != nn:
                    syn_idx.setdefault(ns, []).append((db_id, cname))
                    flat_list.append((ns, db_id, cname, False))

    return name_idx, syn_idx, flat_list


def current_max_node_id() -> int:
    max_id = 0
    for path in (NODE2ID_JSON, APPENDED_JSON):
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        for v in d.values():
            if isinstance(v, int):
                max_id = max(max_id, v)
            elif isinstance(v, dict) and "knowddi_id" in v:
                max_id = max(max_id, v["knowddi_id"])
    return max_id


# ── matching tiers ────────────────────────────────────────────────────────────

def sim(a: str, b: str) -> float:
    return SequenceMatcher(None, a, b).ratio()


def find_drug(
    query: str,
    name_idx: dict,
    syn_idx:  dict,
    flat_list: list,
) -> list[tuple[float, str, str, str]]:
    """
    Returns list of (score, db_id, canonical_name, match_tier) sorted best-first.
    score=1.0 for exact; lower for fuzzy.
    Empty list = nothing found above threshold.
    """
    nq     = norm(query)
    qtoks  = tokenize(query)
    seen: dict[str, tuple[float, str, str]] = {}  # db_id → (score, cname, tier)

    def add(db_id, cname, score, tier):
        if db_id not in seen or score > seen[db_id][0]:
            seen[db_id] = (score, cname, tier)

    # ── tier 1: exact name ───────────────────────────────────────────────────
    for db_id, cname in name_idx.get(nq, []):
        add(db_id, cname, 1.0, "exact_name")

    # ── tier 2: exact synonym ────────────────────────────────────────────────
    for db_id, cname in syn_idx.get(nq, []):
        add(db_id, cname, 1.0, "exact_synonym")

    if seen:
        return _sorted(seen)

    MIN_SUB = 5   # both sides must be at least this long to avoid short-abbrev false matches

    def substr_score(nq: str, nn: str, base: float) -> float:
        """
        Score a substring match.
        - prefix match (candidate starts with query): near-exact, highest score.
          e.g. 'atracurium' is prefix of 'atracurium besylate' → wins over 'cisatracurium'.
        - query buried inside candidate: score by coverage fraction.
        - candidate buried inside query (short word inside long query): heavily penalised.
          e.g. 'calcium' inside 'calcium folinate' → low score, real answer may be elsewhere.
        """
        if nq in nn:
            if nn.startswith(nq):           # query is a prefix of candidate → near-exact
                return base * 0.93
            coverage = len(nq) / len(nn)
            return base * (0.50 + 0.40 * coverage)
        else:  # nn in nq — candidate is a substring of a longer query
            coverage = len(nn) / len(nq)
            return base * 0.50 * coverage   # heavy penalty: generic word inside longer phrase

    # ── tier 3: substring match against canonical names only ─────────────────
    for nn, db_id, cname, is_name in flat_list:
        if not is_name:
            continue
        if len(nq) >= MIN_SUB and len(nn) >= MIN_SUB and (nq in nn or nn in nq):
            add(db_id, cname, substr_score(nq, nn, 0.95), "substr_name")

    # ── tier 4: substring match against synonyms ─────────────────────────────
    for nn, db_id, cname, is_name in flat_list:
        if is_name:
            continue
        if len(nq) >= MIN_SUB and len(nn) >= MIN_SUB and (nq in nn or nn in nq):
            add(db_id, cname, substr_score(nq, nn, 0.90), "substr_synonym")

    if seen:
        return _sorted(seen)

    # ── tiers 5+6: token overlap + fuzzy (scan full flat list) ───────────────
    for nn, db_id, cname, is_name in flat_list:
        # token overlap
        if qtoks:
            ctoks     = set(nn.split())
            overlap   = len(qtoks & ctoks) / len(qtoks)
            if overlap >= TOKEN_MIN_OVERLAP:
                tier  = "token_name" if is_name else "token_synonym"
                score = 0.80 + 0.10 * overlap
                add(db_id, cname, score, tier)

        # fuzzy character similarity
        ratio = sim(nq, nn)
        if ratio >= FUZZY_THRESHOLD:
            tier  = "fuzzy_name" if is_name else "fuzzy_synonym"
            score = ratio * (0.95 if is_name else 0.90)
            add(db_id, cname, score, tier)

    return _sorted(seen)



def _sorted(seen: dict) -> list:
    return sorted(
        [(score, db_id, cname, tier) for db_id, (score, cname, tier) in seen.items()],
        key=lambda x: -x[0]
    )


# ── pretty tier label ─────────────────────────────────────────────────────────

HIGH_CONFIDENCE = {"exact_name", "exact_synonym", "substr_name", "substr_synonym"}

def is_high_confidence(tier: str) -> bool:
    return tier in HIGH_CONFIDENCE


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("drugs", nargs="*", help="Drug names as direct arguments")
    parser.add_argument("--file", "-f",
                        help="Input file: plain text (one name per line) OR a JSON list/dict of drug names")
    parser.add_argument("--skip", "-s", nargs="+", metavar="JSON_FILE",
                        help="JSON mapping files whose keys to skip (e.g. MAPPED_286.json APPENDED_7.json). "
                             "Drugs already in these files are excluded from the search.")
    parser.add_argument("--interactive", "-I", action="store_true",
                        help="Prompt to confirm fuzzy/token matches interactively")
    parser.add_argument("--top", type=int, default=5,
                        help="Number of fuzzy candidates to show (default 5)")
    args = parser.parse_args()

    drug_names: list[str] = list(args.drugs)

    # ── load from file (plain text or JSON list/dict) ─────────────────────────
    if args.file:
        path = Path(args.file)
        if not path.exists():
            print(f"File not found: {args.file}")
            sys.exit(1)
        if path.suffix.lower() == ".json":
            with open(path) as fh:
                data = json.load(fh)
            if isinstance(data, list):
                drug_names += [str(x).strip() for x in data if x]
            elif isinstance(data, dict):
                drug_names += [k.strip() for k in data.keys() if k]
            else:
                print("JSON file must contain a list or dict at the top level.")
                sys.exit(1)
        else:
            with open(path) as fh:
                drug_names += [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    # ── build set of already-mapped drug names to skip ────────────────────────
    already_mapped: set[str] = set()
    for skip_path in (args.skip or []):
        sp = Path(skip_path)
        if not sp.exists():
            print(f"  ⚠ skip file not found, ignored: {skip_path}")
            continue
        with open(sp) as fh:
            d = json.load(fh)
        if isinstance(d, dict):
            already_mapped.update(k.lower().strip() for k in d.keys())
        elif isinstance(d, list):
            already_mapped.update(str(x).lower().strip() for x in d)

    if already_mapped:
        before = len(drug_names)
        drug_names = [n for n in drug_names if n.lower().strip() not in already_mapped]
        print(f"Skipped {before - len(drug_names)} already-mapped drugs "
              f"({before} → {len(drug_names)} remaining)\n")

    if not drug_names:
        print("No new drug names to process (all already mapped or none provided).")
        sys.exit(0)

    print(f"\nLoading DrugBank index …")
    name_idx, syn_idx, flat_list = load_index()
    next_id = current_max_node_id() + 1
    print(f"Next available node ID: {next_id}\n")

    results:   dict        = {}
    not_found: list[str]   = []
    need_review: list[str] = []

    for drug in drug_names:
        candidates = find_drug(drug.lower(), name_idx, syn_idx, flat_list)
        candidates = candidates[:args.top]

        if not candidates:
            print(f"  NOT FOUND  : {drug}")
            not_found.append(drug)
            continue

        best_score, best_db_id, best_cname, best_tier = candidates[0]
        high = is_high_confidence(best_tier)

        if high:
            # auto-accept
            results[drug.lower()] = {
                "drugbank_id":   best_db_id,
                "knowddi_id":    next_id,
                "match":         best_tier,
                "drugbank_name": best_cname,
            }
            print(f"  MAPPED     : {drug!r:45s} → {best_db_id}  '{best_cname}'  [{best_tier}]  node={next_id}")
            if len(candidates) > 1:
                print(f"    alternatives: " +
                      ", ".join(f"{c[1]} '{c[2]}'" for c in candidates[1:3]))
            next_id += 1

        else:
            # low-confidence: show candidates
            print(f"\n  ⚠ FUZZY    : {drug!r}")
            for i, (sc, db_id, cname, tier) in enumerate(candidates):
                print(f"    [{i+1}] {db_id}  '{cname}'  [{tier}]  score={sc:.2f}")

            if args.interactive:
                choice = input(
                    f"  Pick 1–{len(candidates)}, 0=skip, or type a DrugBank ID (e.g. DB00123): "
                ).strip()
                if choice.startswith("DB") or choice.startswith("db"):
                    # manual entry
                    manual_id = choice.upper()
                    results[drug.lower()] = {
                        "drugbank_id":   manual_id,
                        "knowddi_id":    next_id,
                        "match":         "manual",
                        "drugbank_name": drug,
                    }
                    print(f"    Saved manual: {manual_id}  node={next_id}")
                    next_id += 1
                elif choice.isdigit() and 1 <= int(choice) <= len(candidates):
                    idx = int(choice) - 1
                    _, sel_db_id, sel_cname, sel_tier = candidates[idx]
                    results[drug.lower()] = {
                        "drugbank_id":   sel_db_id,
                        "knowddi_id":    next_id,
                        "match":         sel_tier,
                        "drugbank_name": sel_cname,
                    }
                    print(f"    Confirmed: {sel_db_id}  '{sel_cname}'  node={next_id}")
                    next_id += 1
                else:
                    print(f"    Skipped — add to output manually if needed.")
                    need_review.append(drug)
            else:
                # non-interactive: write best fuzzy match but flag it
                results[drug.lower()] = {
                    "drugbank_id":   best_db_id,
                    "knowddi_id":    next_id,
                    "match":         best_tier,
                    "drugbank_name": best_cname,
                    "NEEDS_REVIEW":  True,
                }
                print(f"    Auto-wrote best fuzzy match [{best_tier}] — verify manually.")
                need_review.append(drug)
                next_id += 1
            print()

    print(f"\n{'─'*60}")
    print(f"Mapped        : {len(results)}")
    print(f"Not found     : {len(not_found)}  → {not_found}")
    if need_review:
        print(f"Needs review  : {len(need_review)}  → {need_review}")
        print("  Re-run with --interactive to confirm fuzzy matches one-by-one.")

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved → {OUT_JSON}")


if __name__ == "__main__":
    main()
