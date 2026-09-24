#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re
import time

import pandas as pd
import requests


# ============================================================
# pfam.py
#
# Pfam membership vrstva pro Lucka-Structures.
#
# Použití:
#   python pfam.py --pfam PF00407
#   python pfam.py project7 --pfam PF00407
#
# Data čte z InterPro REST API a ukládá do:
#   data/processed/projectN/pfam/
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Stáhne všechny UniProtKB proteiny odpovídající zadanému Pfam entry "
        "a připraví taxonomickou evidence vrstvu pro viewer."
    )
)
parser.add_argument(
    "project",
    nargs="?",
    help="Např. project7. Bez argumentu se použije data/current_project.txt."
)
parser.add_argument(
    "--pfam",
    required=True,
    help="Pfam accession, např. PF00407."
)
parser.add_argument(
    "--reviewed-only",
    action="store_true",
    help="Použít jen UniProtKB/Swiss-Prot. Default je reviewed + unreviewed UniProtKB."
)
parser.add_argument(
    "--page-size",
    type=int,
    default=200,
    help="Velikost stránky InterPro API (default 200; maximum bývá 200)."
)
parser.add_argument(
    "--sleep",
    type=float,
    default=0.15,
    help="Pauza mezi stránkami API v sekundách (default 0.15)."
)
args = parser.parse_args()

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
PROCESSED_ROOT_DIR = DATA_DIR / "processed"
CURRENT_PROJECT_FILE = DATA_DIR / "current_project.txt"

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
HEADERS = {
    "Accept": "application/json",
    "User-Agent": "Lucka-Structures Pfam membership layer (research/educational use)",
}
PROJECT_RE = re.compile(r"^project\d+$", re.IGNORECASE)
PFAM_RE = re.compile(r"^PF\d{5}$", re.IGNORECASE)


def get_project_id():
    if args.project:
        project_id = args.project.strip().lower()
    else:
        if not CURRENT_PROJECT_FILE.exists():
            raise FileNotFoundError(
                "Neexistuje data/current_project.txt. Zadej projectN ručně."
            )
        project_id = CURRENT_PROJECT_FILE.read_text(encoding="utf-8").strip().lower()

    if not PROJECT_RE.fullmatch(project_id):
        raise ValueError(f"Neplatné project ID: {project_id!r}")
    return project_id


def request_json(session, url, params=None, attempts=5):
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                url,
                params=params,
                headers=HEADERS,
                timeout=90,
            )

            if response.status_code == 204:
                return {}

            if response.status_code == 429 or response.status_code >= 500:
                wait = min(30, 2 ** attempt)
                print(f"  API HTTP {response.status_code}; opakuji za {wait} s…")
                time.sleep(wait)
                continue

            response.raise_for_status()
            return response.json()
        except Exception as error:
            last_error = error
            if attempt == attempts:
                break
            wait = min(30, 2 ** attempt)
            print(f"  API chyba: {error}; opakuji za {wait} s…")
            time.sleep(wait)

    raise RuntimeError(f"InterPro API selhalo po {attempts} pokusech: {last_error}")


project_id = get_project_id()
pfam_id = args.pfam.strip().upper()
if not PFAM_RE.fullmatch(pfam_id):
    raise ValueError(
        f"Neplatné Pfam accession {pfam_id!r}. Očekávám např. PF00407."
    )

processed_dir = PROCESSED_ROOT_DIR / project_id
project_json_file = processed_dir / "project.json"

if not processed_dir.exists():
    raise FileNotFoundError(f"Projekt neexistuje: {processed_dir}")

pfam_dir = processed_dir / "pfam"
pfam_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print(f"PROJECT: {project_id}")
print(f"PFAM:    {pfam_id}")
print(f"MODE:    {'reviewed only' if args.reviewed_only else 'all UniProtKB'}")
print("=" * 70)

session = requests.Session()

# ------------------------------------------------------------
# Metadata Pfam entry
# ------------------------------------------------------------
entry_url = f"{INTERPRO_API}/entry/pfam/{pfam_id}/"
entry_payload = request_json(session, entry_url)
entry_meta = entry_payload.get("metadata", {}) if isinstance(entry_payload, dict) else {}
name_field = entry_meta.get("name")
if isinstance(name_field, dict):
    pfam_name = name_field.get("name") or name_field.get("short") or pfam_id
else:
    pfam_name = name_field or pfam_id

integrated_interpro = entry_meta.get("integrated")
print(f"Pfam název: {pfam_name}")
if integrated_interpro:
    print(f"Integrované InterPro: {integrated_interpro}")

# ------------------------------------------------------------
# Všechny proteiny odpovídající entry
# ------------------------------------------------------------
protein_source = "reviewed" if args.reviewed_only else "uniprot"
url = f"{INTERPRO_API}/protein/{protein_source}/entry/pfam/{pfam_id}/"
params = {"page_size": max(1, min(args.page_size, 200))}

rows = []
seen_accessions = set()
page = 0
expected_count = None

while url:
    page += 1
    payload = request_json(session, url, params=params if page == 1 else None)
    params = None

    if expected_count is None:
        expected_count = payload.get("count")
        if expected_count is not None:
            print(f"InterPro hlásí {expected_count} proteinů.")

    results = payload.get("results", []) or []

    for item in results:
        meta = item.get("metadata", {}) or {}
        accession = str(meta.get("accession") or "").strip()
        if not accession or accession in seen_accessions:
            continue
        seen_accessions.add(accession)

        organism = meta.get("source_organism") or {}
        taxid_raw = organism.get("taxId")
        try:
            taxid = int(taxid_raw) if taxid_raw not in (None, "") else None
        except (TypeError, ValueError):
            taxid = None

        gene = meta.get("gene")
        if isinstance(gene, dict):
            gene = gene.get("name") or gene.get("id") or json.dumps(gene, ensure_ascii=False)
        elif isinstance(gene, list):
            gene = ", ".join(str(x) for x in gene)

        rows.append({
            "uniprot": accession,
            "protein_name": meta.get("name"),
            "organism": organism.get("scientificName") or organism.get("fullName"),
            "taxid": taxid,
            "gene": gene,
            "length": meta.get("length"),
            "source_database": meta.get("source_database"),
            "pfam_id": pfam_id,
            "pfam_name": pfam_name,
        })

    print(
        f"  stránka {page}: +{len(results)} záznamů; "
        f"unikátních proteinů {len(rows)}"
    )

    url = payload.get("next")
    if url:
        time.sleep(max(0.0, args.sleep))

members_df = pd.DataFrame(rows)
if members_df.empty:
    members_df = pd.DataFrame(columns=[
        "uniprot", "protein_name", "organism", "taxid", "gene", "length",
        "source_database", "pfam_id", "pfam_name"
    ])

members_file = pfam_dir / "pfam_members.tsv"
members_df.to_csv(members_file, sep="\t", index=False)

with_taxid = members_df[members_df["taxid"].notna()].copy()
if not with_taxid.empty:
    with_taxid["taxid"] = with_taxid["taxid"].astype(int)

    taxa_rows = []
    for taxid, group in with_taxid.groupby("taxid", sort=False):
        organisms = [x for x in group["organism"].dropna().astype(str).unique() if x]
        accessions = group["uniprot"].dropna().astype(str).tolist()
        taxa_rows.append({
            "taxid": int(taxid),
            "organism": organisms[0] if organisms else "",
            "family_member_count": int(len(group)),
            "pfam_id": pfam_id,
            "pfam_name": pfam_name,
            "uniprot_accessions": ",".join(accessions),
        })

    taxa_df = pd.DataFrame(taxa_rows).sort_values(
        ["family_member_count", "organism"], ascending=[False, True]
    )
else:
    taxa_df = pd.DataFrame(columns=[
        "taxid", "organism", "family_member_count", "pfam_id", "pfam_name",
        "uniprot_accessions"
    ])

taxa_file = pfam_dir / "pfam_taxa.tsv"
taxa_df.to_csv(taxa_file, sep="\t", index=False)

# Viewer JSON
viewer_taxa = {}
for _, row in taxa_df.iterrows():
    taxid = int(row["taxid"])
    accessions = [
        x for x in str(row.get("uniprot_accessions") or "").split(",") if x
    ]
    viewer_taxa[str(taxid)] = {
        "taxid": taxid,
        "organism": row.get("organism") or "",
        "family_member": True,
        "family_member_count": int(row.get("family_member_count") or 0),
        "pfam_id": pfam_id,
        "pfam_name": pfam_name,
        "uniprot_accessions": accessions,
    }

evidence = {
    "source": "InterPro API / Pfam",
    "pfam_id": pfam_id,
    "pfam_name": pfam_name,
    "integrated_interpro": integrated_interpro,
    "reviewed_only": bool(args.reviewed_only),
    "taxa": viewer_taxa,
}
evidence_file = pfam_dir / "pfam_evidence.json"
evidence_file.write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")

summary = {
    "project_id": project_id,
    "pfam_id": pfam_id,
    "pfam_name": pfam_name,
    "integrated_interpro": integrated_interpro,
    "reviewed_only": bool(args.reviewed_only),
    "api_reported_count": expected_count,
    "members": {
        "count": int(len(members_df)),
        "with_taxid": int(members_df["taxid"].notna().sum()) if "taxid" in members_df else 0,
        "without_taxid": int(members_df["taxid"].isna().sum()) if "taxid" in members_df else 0,
        "unique_taxids": int(taxa_df["taxid"].nunique()) if not taxa_df.empty else 0,
    },
    "files": {
        "members_tsv": "pfam/pfam_members.tsv",
        "taxa_tsv": "pfam/pfam_taxa.tsv",
        "evidence_json": "pfam/pfam_evidence.json",
    },
}
summary_file = pfam_dir / "pfam_summary.json"
summary_file.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

# project.json metadata
if project_json_file.exists():
    project_meta = json.loads(project_json_file.read_text(encoding="utf-8"))
else:
    project_meta = {"project_id": project_id}

project_meta["pfam"] = {
    "pfam_id": pfam_id,
    "pfam_name": pfam_name,
    "integrated_interpro": integrated_interpro,
    "reviewed_only": bool(args.reviewed_only),
    "summary": "pfam/pfam_summary.json",
    "members_tsv": "pfam/pfam_members.tsv",
    "taxa_tsv": "pfam/pfam_taxa.tsv",
    "evidence_json": "pfam/pfam_evidence.json",
}
project_json_file.write_text(
    json.dumps(project_meta, ensure_ascii=False, indent=2), encoding="utf-8"
)

print()
print("Pfam evidence hotová:")
print(f"  proteinů: {len(members_df)}")
print(f"  s Tax ID: {members_df['taxid'].notna().sum() if 'taxid' in members_df else 0}")
print(f"  taxonů:   {len(taxa_df)}")
print(f"  {members_file}")
print(f"  {taxa_file}")
print(f"  {evidence_file}")
