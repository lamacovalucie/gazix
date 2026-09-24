#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re
import time

import pandas as pd
import requests


# ============================================================
# family.py
#
# Obecná InterPro family-membership vrstva pro Lucka-Structures.
#
# Použití:
#   python family.py --interpro IPR003571
#   python family.py project7 --interpro IPR003571
#
# Výstupy:
#   data/processed/projectN/family/
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Stáhne všechny UniProtKB proteiny patřící do zadaného "
        "InterPro entry a připraví taxonomickou evidence vrstvu."
    )
)
parser.add_argument(
    "project",
    nargs="?",
    help="Např. project7. Bez argumentu se použije data/current_project.txt."
)
parser.add_argument(
    "--interpro",
    required=True,
    help="InterPro accession, např. IPR003571."
)
parser.add_argument(
    "--reviewed-only",
    action="store_true",
    help="Použít jen reviewed UniProtKB/Swiss-Prot."
)
parser.add_argument(
    "--page-size",
    type=int,
    default=200,
    help="Velikost stránky InterPro API (default 200)."
)
parser.add_argument(
    "--sleep",
    type=float,
    default=0.15,
    help="Pauza mezi stránkami API v sekundách."
)
parser.add_argument(
    "--max-attempts",
    type=int,
    default=12,
    help=(
        "Kolikrát opakovat dočasně neúspěšný InterPro request "
        "(default 12)."
    )
)
parser.add_argument(
    "--fresh",
    action="store_true",
    help=(
        "Ignorovat rozpracovaný checkpoint a stáhnout InterPro "
        "evidence od začátku."
    )
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
    "User-Agent": "Lucka-Structures InterPro family layer (research/educational use)",
}

PROJECT_RE = re.compile(r"^project\d+$", re.IGNORECASE)
INTERPRO_RE = re.compile(r"^IPR\d{6}$", re.IGNORECASE)


def get_project_id():
    if args.project:
        project_id = args.project.strip().lower()
    else:
        if not CURRENT_PROJECT_FILE.exists():
            raise FileNotFoundError(
                "Neexistuje data/current_project.txt. Zadej projectN ručně."
            )
        project_id = CURRENT_PROJECT_FILE.read_text(
            encoding="utf-8"
        ).strip().lower()

    if not PROJECT_RE.fullmatch(project_id):
        raise ValueError(f"Neplatné project ID: {project_id!r}")

    return project_id


def request_json(session, url, params=None, attempts=None):
    """GET JSON s odolným retry pro přechodné chyby InterPro API."""

    if attempts is None:
        attempts = max(1, int(args.max_attempts))

    last_error = None

    # 408 je u velkých InterPro dotazů běžný timeout serveru.
    retryable_statuses = {
        408, 425, 429,
        500, 502, 503, 504,
        520, 521, 522, 523, 524,
    }

    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                url,
                params=params,
                headers=HEADERS,
                timeout=120,
            )

            if response.status_code == 204:
                return {}

            if response.status_code in retryable_statuses:
                retry_after = response.headers.get("Retry-After")

                try:
                    wait = float(retry_after) if retry_after else None
                except (TypeError, ValueError):
                    wait = None

                if wait is None:
                    # 2, 4, 8, 16, 32, 60, 60, ...
                    wait = min(60, 2 ** attempt)

                last_error = requests.HTTPError(
                    f"{response.status_code} {response.reason or ''}".strip()
                )

                if attempt == attempts:
                    break

                print(
                    f"  API HTTP {response.status_code}; "
                    f"pokus {attempt}/{attempts}, "
                    f"opakuji za {wait:g} s…"
                )
                time.sleep(wait)
                continue

            response.raise_for_status()

            try:
                return response.json()
            except ValueError as error:
                preview = (response.text or "")[:200].replace("\\n", " ")
                last_error = error

                if attempt == attempts:
                    break

                wait = min(60, 2 ** attempt)
                print(
                    "  API nevrátilo platný JSON "
                    f"(odpověď: {preview!r}); "
                    f"pokus {attempt}/{attempts}, "
                    f"opakuji za {wait:g} s…"
                )
                time.sleep(wait)
                continue

        except requests.RequestException as error:
            last_error = error

            if attempt == attempts:
                break

            wait = min(60, 2 ** attempt)
            print(
                f"  API chyba: {error}; "
                f"pokus {attempt}/{attempts}, "
                f"opakuji za {wait:g} s…"
            )
            time.sleep(wait)

    raise RuntimeError(
        f"InterPro API selhalo po {attempts} pokusech: {last_error}"
    )


project_id = get_project_id()
interpro_id = args.interpro.strip().upper()

if not INTERPRO_RE.fullmatch(interpro_id):
    raise ValueError(
        f"Neplatné InterPro accession {interpro_id!r}. "
        "Očekávám např. IPR003571."
    )

processed_dir = PROCESSED_ROOT_DIR / project_id
project_json_file = processed_dir / "project.json"

if not processed_dir.exists():
    raise FileNotFoundError(f"Projekt neexistuje: {processed_dir}")

family_dir = processed_dir / "family"
family_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print(f"PROJECT:  {project_id}")
print(f"INTERPRO: {interpro_id}")
print(
    f"MODE:     "
    f"{'reviewed only' if args.reviewed_only else 'all UniProtKB'}"
)
print("=" * 70)

session = requests.Session()

# ------------------------------------------------------------
# Metadata InterPro entry
# ------------------------------------------------------------

entry_url = f"{INTERPRO_API}/entry/interpro/{interpro_id}/"
entry_payload = request_json(session, entry_url)

entry_meta = (
    entry_payload.get("metadata", {})
    if isinstance(entry_payload, dict)
    else {}
)

name_field = entry_meta.get("name")

if isinstance(name_field, dict):
    family_name = (
        name_field.get("name")
        or name_field.get("short")
        or interpro_id
    )
else:
    family_name = name_field or interpro_id

entry_type = entry_meta.get("type") or entry_meta.get("entry_type")

print(f"InterPro název: {family_name}")
if entry_type:
    print(f"Typ entry: {entry_type}")

# ------------------------------------------------------------
# Proteiny patřící do InterPro entry
# ------------------------------------------------------------

protein_source = "reviewed" if args.reviewed_only else "uniprot"

initial_url = (
    f"{INTERPRO_API}/protein/{protein_source}/"
    f"entry/interpro/{interpro_id}/"
)

page_size = max(1, min(args.page_size, 200))
params = {
    "page_size": page_size
}

# Checkpoint soubory. Partial TSV se průběžně přepisuje atomicky
# po každé úspěšně stažené stránce.
checkpoint_file = family_dir / "family_download_checkpoint.json"
partial_members_file = family_dir / "family_members.partial.tsv"

rows = []
seen_accessions = set()
page = 0
expected_count = None
url = initial_url
resume_from_checkpoint = False


def save_partial_members(rows_to_save):
    partial_df = pd.DataFrame(rows_to_save)

    if partial_df.empty:
        partial_df = pd.DataFrame(
            columns=[
                "uniprot",
                "protein_name",
                "organism",
                "taxid",
                "gene",
                "length",
                "source_database",
                "interpro_id",
                "interpro_name",
            ]
        )

    tmp_file = partial_members_file.with_suffix(
        partial_members_file.suffix + ".tmp"
    )

    partial_df.to_csv(
        tmp_file,
        sep="\t",
        index=False
    )

    tmp_file.replace(partial_members_file)


def save_checkpoint(next_url, completed_page):
    checkpoint = {
        "project_id": project_id,
        "interpro_id": interpro_id,
        "reviewed_only": bool(args.reviewed_only),
        "page_size": page_size,
        "completed_page": int(completed_page),
        "expected_count": expected_count,
        "next_url": next_url,
        "partial_members_file": partial_members_file.name,
        "saved_members": len(rows),
    }

    tmp_file = checkpoint_file.with_suffix(
        checkpoint_file.suffix + ".tmp"
    )

    tmp_file.write_text(
        json.dumps(
            checkpoint,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    tmp_file.replace(checkpoint_file)


# ------------------------------------------------------------
# Případné pokračování z checkpointu
# ------------------------------------------------------------

if args.fresh:
    for stale_file in (
        checkpoint_file,
        partial_members_file,
    ):
        if stale_file.exists():
            stale_file.unlink()

elif checkpoint_file.exists() and partial_members_file.exists():
    try:
        checkpoint = json.loads(
            checkpoint_file.read_text(
                encoding="utf-8"
            )
        )

        checkpoint_matches = (
            checkpoint.get("project_id") == project_id
            and checkpoint.get("interpro_id") == interpro_id
            and bool(checkpoint.get("reviewed_only"))
                == bool(args.reviewed_only)
            and int(checkpoint.get("page_size", -1))
                == page_size
            and checkpoint.get("next_url")
        )

        if checkpoint_matches:
            partial_df = pd.read_csv(
                partial_members_file,
                sep="\t"
            )

            rows = partial_df.where(
                pd.notnull(partial_df),
                None
            ).to_dict(
                orient="records"
            )

            seen_accessions = {
                str(row.get("uniprot") or "").strip()
                for row in rows
                if str(row.get("uniprot") or "").strip()
            }

            page = int(
                checkpoint.get(
                    "completed_page",
                    0
                )
            )

            expected_count = checkpoint.get(
                "expected_count"
            )

            url = checkpoint.get(
                "next_url"
            )

            # Cursor URL už page_size obsahuje, proto při resume
            # neposíláme params znovu.
            params = None
            resume_from_checkpoint = True

            print()
            print(
                "Nalezen rozpracovaný InterPro download — "
                "pokračuji z checkpointu."
            )
            print(
                f"  hotových stránek: {page}"
            )
            print(
                f"  uložených proteinů: {len(rows)}"
            )
            if expected_count is not None:
                print(
                    f"  InterPro celkem: {expected_count}"
                )
            print()

        else:
            print(
                "Checkpoint neodpovídá aktuálnímu projektu/InterPro "
                "nastavení; stahuji od začátku."
            )

    except Exception as error:
        print(
            "Checkpoint se nepodařilo načíst; "
            f"stahuji od začátku ({error})."
        )
        rows = []
        seen_accessions = set()
        page = 0
        expected_count = None
        url = initial_url
        params = {
            "page_size": page_size
        }


while url:
    page += 1

    try:
        payload = request_json(
            session,
            url,
            params=params if page == 1 and not resume_from_checkpoint else None
        )
    except Exception:
        # Před vyhozením chyby znovu uložíme poslední bezpečný stav.
        # next_url zde stále ukazuje na stránku, která právě selhala.
        save_partial_members(rows)
        save_checkpoint(url, page - 1)

        print()
        print(
            "Download byl přerušen, ale dosavadní stav je uložen."
        )
        print(
            f"Po dalším spuštění stejného příkazu budu pokračovat "
            f"od stránky {page}."
        )
        print(
            f"Checkpoint: {checkpoint_file}"
        )
        raise

    # params se používají pouze u úplně první stránky.
    params = None
    resume_from_checkpoint = False

    if expected_count is None:
        expected_count = payload.get("count")
        if expected_count is not None:
            print(f"InterPro hlásí {expected_count} proteinů.")

    results = payload.get("results", []) or []

    for item in results:
        meta = item.get("metadata", {}) or {}

        accession = str(
            meta.get("accession") or ""
        ).strip()

        if not accession or accession in seen_accessions:
            continue

        seen_accessions.add(accession)

        organism = meta.get("source_organism") or {}

        taxid_raw = organism.get("taxId")

        try:
            taxid = (
                int(taxid_raw)
                if taxid_raw not in (None, "")
                else None
            )
        except (TypeError, ValueError):
            taxid = None

        gene = meta.get("gene")

        if isinstance(gene, dict):
            gene = (
                gene.get("name")
                or gene.get("id")
                or json.dumps(
                    gene,
                    ensure_ascii=False
                )
            )

        elif isinstance(gene, list):
            gene = ", ".join(
                str(x)
                for x in gene
            )

        rows.append({
            "uniprot": accession,
            "protein_name": meta.get("name"),
            "organism": (
                organism.get("scientificName")
                or organism.get("fullName")
            ),
            "taxid": taxid,
            "gene": gene,
            "length": meta.get("length"),
            "source_database": meta.get("source_database"),
            "interpro_id": interpro_id,
            "interpro_name": family_name,
        })

    next_url = payload.get("next")

    # Každá úspěšná stránka se uloží. Když nás EBI později shodí,
    # příští běh začne přesně od next_url.
    save_partial_members(rows)
    save_checkpoint(next_url, page)

    print(
        f"  stránka {page}: +{len(results)} záznamů; "
        f"unikátních proteinů {len(rows)} "
        f"[checkpoint uložen]"
    )

    url = next_url

    if url:
        time.sleep(
            max(0.0, args.sleep)
        )


members_df = pd.DataFrame(rows)

if members_df.empty:
    members_df = pd.DataFrame(
        columns=[
            "uniprot",
            "protein_name",
            "organism",
            "taxid",
            "gene",
            "length",
            "source_database",
            "interpro_id",
            "interpro_name",
        ]
    )

members_file = (
    family_dir / "family_members.tsv"
)

members_df.to_csv(
    members_file,
    sep="\t",
    index=False
)

# Finální members TSV je bezpečně zapsaný; rozpracovaný checkpoint už
# není potřeba. Při pádu před tímto místem zůstane zachován.
if checkpoint_file.exists():
    checkpoint_file.unlink()

if partial_members_file.exists():
    partial_members_file.unlink()

with_taxid = members_df[
    members_df["taxid"].notna()
].copy()

if not with_taxid.empty:

    with_taxid["taxid"] = (
        with_taxid["taxid"]
        .astype(int)
    )

    taxa_rows = []

    for taxid, group in with_taxid.groupby(
        "taxid",
        sort=False
    ):

        organisms = [
            x
            for x in (
                group["organism"]
                .dropna()
                .astype(str)
                .unique()
            )
            if x
        ]

        accessions = (
            group["uniprot"]
            .dropna()
            .astype(str)
            .tolist()
        )

        taxa_rows.append({
            "taxid": int(taxid),
            "organism": (
                organisms[0]
                if organisms
                else ""
            ),
            "family_member_count": int(
                len(group)
            ),
            "interpro_id": interpro_id,
            "interpro_name": family_name,
            "uniprot_accessions": ",".join(
                accessions
            ),
        })

    taxa_df = pd.DataFrame(
        taxa_rows
    ).sort_values(
        [
            "family_member_count",
            "organism"
        ],
        ascending=[
            False,
            True
        ]
    )

else:

    taxa_df = pd.DataFrame(
        columns=[
            "taxid",
            "organism",
            "family_member_count",
            "interpro_id",
            "interpro_name",
            "uniprot_accessions",
        ]
    )


taxa_file = (
    family_dir / "family_taxa.tsv"
)

taxa_df.to_csv(
    taxa_file,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Viewer JSON
# ------------------------------------------------------------

viewer_taxa = {}

for _, row in taxa_df.iterrows():

    taxid = int(
        row["taxid"]
    )

    accessions = [
        x
        for x in str(
            row.get(
                "uniprot_accessions"
            ) or ""
        ).split(",")
        if x
    ]

    viewer_taxa[
        str(taxid)
    ] = {
        "taxid": taxid,
        "organism": (
            row.get("organism")
            or ""
        ),
        "family_member": True,
        "family_member_count": int(
            row.get(
                "family_member_count"
            ) or 0
        ),
        "interpro_id": interpro_id,
        "interpro_name": family_name,
        "uniprot_accessions": accessions,
    }


evidence = {
    "source": "InterPro API",
    "interpro_id": interpro_id,
    "interpro_name": family_name,
    "entry_type": entry_type,
    "reviewed_only": bool(
        args.reviewed_only
    ),
    "taxa": viewer_taxa,
}

evidence_file = (
    family_dir
    / "family_evidence.json"
)

evidence_file.write_text(
    json.dumps(
        evidence,
        ensure_ascii=False,
        indent=2
    ),
    encoding="utf-8"
)


summary = {
    "project_id": project_id,
    "interpro_id": interpro_id,
    "interpro_name": family_name,
    "entry_type": entry_type,
    "reviewed_only": bool(
        args.reviewed_only
    ),
    "api_reported_count": expected_count,
    "members": {
        "count": int(
            len(members_df)
        ),
        "with_taxid": int(
            members_df[
                "taxid"
            ].notna().sum()
        ) if "taxid" in members_df else 0,
        "without_taxid": int(
            members_df[
                "taxid"
            ].isna().sum()
        ) if "taxid" in members_df else 0,
        "unique_taxids": int(
            taxa_df[
                "taxid"
            ].nunique()
        ) if not taxa_df.empty else 0,
    },
    "files": {
        "members_tsv": (
            "family/family_members.tsv"
        ),
        "taxa_tsv": (
            "family/family_taxa.tsv"
        ),
        "evidence_json": (
            "family/family_evidence.json"
        ),
    },
}

summary_file = (
    family_dir
    / "family_summary.json"
)

summary_file.write_text(
    json.dumps(
        summary,
        ensure_ascii=False,
        indent=2
    ),
    encoding="utf-8"
)

# ------------------------------------------------------------
# project.json metadata
# ------------------------------------------------------------

if project_json_file.exists():

    project_meta = json.loads(
        project_json_file.read_text(
            encoding="utf-8"
        )
    )

else:

    project_meta = {
        "project_id": project_id
    }


project_meta["family"] = {
    "source": "InterPro",
    "interpro_id": interpro_id,
    "interpro_name": family_name,
    "entry_type": entry_type,
    "reviewed_only": bool(
        args.reviewed_only
    ),
    "summary": (
        "family/family_summary.json"
    ),
    "members_tsv": (
        "family/family_members.tsv"
    ),
    "taxa_tsv": (
        "family/family_taxa.tsv"
    ),
    "evidence_json": (
        "family/family_evidence.json"
    ),
}

project_json_file.write_text(
    json.dumps(
        project_meta,
        ensure_ascii=False,
        indent=2
    ),
    encoding="utf-8"
)


print()
print("InterPro family evidence hotová:")
print(
    f"  proteinů: {len(members_df)}"
)
print(
    "  s Tax ID: "
    f"{members_df['taxid'].notna().sum() if 'taxid' in members_df else 0}"
)
print(
    f"  taxonů:   {len(taxa_df)}"
)
print(f"  {members_file}")
print(f"  {taxa_file}")
print(f"  {evidence_file}")
