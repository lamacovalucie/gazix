#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re
import time

import pandas as pd
import requests


# ============================================================
# family.py — InterPro TAXON-FIRST evidence layer
# ============================================================
#
# Místo stahování všech proteinů dané InterPro entry používá přímo
# taxonomický endpoint InterPro API:
#
#   /api/taxonomy/uniprot/entry/interpro/IPRxxxxxx/
#
# Tím získává organismy/taxony, které mají alespoň jeden protein
# odpovídající dané InterPro entry, bez nutnosti stáhnout desetitisíce
# jednotlivých UniProt záznamů.
#
# Použití:
#
#   python family.py --interpro IPR000916
#   python family.py --interpro IPR003571
#   python family.py --interpro IPR045860
#
# Výstupy:
#
#   data/processed/projectN/family/family_taxa.tsv
#   data/processed/projectN/family/family_evidence.json
#   data/processed/projectN/family/family_summary.json
#
# Skript podporuje checkpoint. Po pádu stejné spuštění naváže.
# ============================================================


parser = argparse.ArgumentParser(
    description=(
        "Stáhne TAXONOMICKOU distribuci zadané InterPro entry "
        "a připraví evidence vrstvu pro viewer."
    )
)

parser.add_argument(
    "project",
    nargs="?",
    help=(
        "Např. project17. Bez argumentu se použije "
        "data/current_project.txt."
    )
)

parser.add_argument(
    "--interpro",
    required=True,
    help="InterPro accession, např. IPR003571."
)

parser.add_argument(
    "--page-size",
    type=int,
    default=200,
    help="Velikost stránky InterPro API (default 200, maximum 200)."
)

parser.add_argument(
    "--sleep",
    type=float,
    default=0.2,
    help="Pauza mezi stránkami API v sekundách."
)

parser.add_argument(
    "--max-attempts",
    type=int,
    default=12,
    help="Počet retry při dočasné chybě API (default 12)."
)

parser.add_argument(
    "--fresh",
    action="store_true",
    help="Ignorovat starý checkpoint a stáhnout taxony od začátku."
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
    "User-Agent": (
        "Lucka-Structures InterPro taxonomy layer "
        "(research/educational use)"
    ),
}

PROJECT_RE = re.compile(r"^project\d+$", re.IGNORECASE)
INTERPRO_RE = re.compile(r"^IPR\d{6}$", re.IGNORECASE)


def get_project_id():
    if args.project:
        project_id = args.project.strip().lower()
    else:
        if not CURRENT_PROJECT_FILE.exists():
            raise FileNotFoundError(
                "Neexistuje data/current_project.txt. "
                "Zadej projectN ručně."
            )

        project_id = CURRENT_PROJECT_FILE.read_text(
            encoding="utf-8"
        ).strip().lower()

    if not PROJECT_RE.fullmatch(project_id):
        raise ValueError(
            f"Neplatné project ID: {project_id!r}"
        )

    return project_id


def request_json(session, url, params=None):
    """Odolný GET JSON pro InterPro API."""

    attempts = max(1, int(args.max_attempts))

    retryable_statuses = {
        408, 425, 429,
        500, 502, 503, 504,
        520, 521, 522, 523, 524,
    }

    last_error = None

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
                last_error = requests.HTTPError(
                    f"{response.status_code} "
                    f"{response.reason or ''}".strip()
                )

                if attempt == attempts:
                    break

                retry_after = response.headers.get(
                    "Retry-After"
                )

                try:
                    wait = (
                        float(retry_after)
                        if retry_after
                        else None
                    )
                except (TypeError, ValueError):
                    wait = None

                if wait is None:
                    wait = min(60, 2 ** attempt)

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
                last_error = error

                if attempt == attempts:
                    break

                wait = min(60, 2 ** attempt)

                preview = (
                    response.text or ""
                )[:200].replace("\n", " ")

                print(
                    "  API nevrátilo platný JSON "
                    f"({preview!r}); "
                    f"opakuji za {wait:g} s…"
                )

                time.sleep(wait)

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
        f"InterPro API selhalo po {attempts} pokusech: "
        f"{last_error}"
    )


def first_nonempty(*values):
    for value in values:
        if value not in (None, "", [], {}):
            return value
    return None


def parse_taxid(item):
    """Tax ID z různých variant InterPro taxonomy JSON."""

    meta = item.get("metadata", {}) or {}

    candidates = [
        meta.get("accession"),
        meta.get("tax_id"),
        meta.get("taxId"),
        meta.get("id"),
        item.get("accession"),
        item.get("tax_id"),
        item.get("taxId"),
        item.get("id"),
    ]

    for value in candidates:
        try:
            if value not in (None, ""):
                return int(value)
        except (TypeError, ValueError):
            continue

    return None


def parse_scientific_name(item):
    meta = item.get("metadata", {}) or {}

    scientific = first_nonempty(
        meta.get("scientific_name"),
        meta.get("scientificName"),
        item.get("scientific_name"),
        item.get("scientificName"),
    )

    if scientific:
        return str(scientific)

    full_name = first_nonempty(
        meta.get("full_name"),
        meta.get("fullName"),
        item.get("full_name"),
        item.get("fullName"),
    )

    if isinstance(full_name, str):
        return full_name

    if isinstance(full_name, dict):
        return str(
            first_nonempty(
                full_name.get("name"),
                full_name.get("scientific"),
                full_name.get("scientific_name"),
                ""
            )
        )

    return ""


def parse_rank(item):
    meta = item.get("metadata", {}) or {}

    return str(
        first_nonempty(
            meta.get("rank"),
            item.get("rank"),
            ""
        )
    )


def find_protein_count(item):
    """
    Pokusí se vytáhnout počet proteinů pro taxon, pokud ho taxonomy
    endpoint poskytuje. Není-li dostupný, vrátí None.
    """

    meta = item.get("metadata", {}) or {}

    direct_candidates = [
        item.get("protein_count"),
        item.get("proteins"),
        meta.get("protein_count"),
        meta.get("proteins"),
    ]

    for value in direct_candidates:
        if isinstance(value, (int, float)):
            return int(value)

    # InterPro často vrací agregované counters v různých tvarech.
    containers = [
        item.get("counters"),
        meta.get("counters"),
    ]

    def walk(value):
        if isinstance(value, dict):
            # Nejdřív explicitní klíče.
            for key in (
                "proteins",
                "protein",
                "protein_count",
                "count",
            ):
                candidate = value.get(key)
                if isinstance(candidate, (int, float)):
                    return int(candidate)

            for child in value.values():
                found = walk(child)
                if found is not None:
                    return found

        elif isinstance(value, list):
            for child in value:
                found = walk(child)
                if found is not None:
                    return found

        return None

    for container in containers:
        found = walk(container)
        if found is not None:
            return found

    return None


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
    raise FileNotFoundError(
        f"Projekt neexistuje: {processed_dir}"
    )

family_dir = processed_dir / "family"
family_dir.mkdir(
    parents=True,
    exist_ok=True
)

page_size = max(
    1,
    min(int(args.page_size), 200)
)

checkpoint_file = (
    family_dir
    / "family_taxonomy_checkpoint.json"
)

partial_file = (
    family_dir
    / "family_taxa.partial.tsv"
)

taxa_file = (
    family_dir
    / "family_taxa.tsv"
)

evidence_file = (
    family_dir
    / "family_evidence.json"
)

summary_file = (
    family_dir
    / "family_summary.json"
)


print("=" * 70)
print(f"PROJECT:  {project_id}")
print(f"INTERPRO: {interpro_id}")
print("MODE:     taxonomy-first")
print("=" * 70)


session = requests.Session()


# ------------------------------------------------------------
# Metadata InterPro entry
# ------------------------------------------------------------

entry_url = (
    f"{INTERPRO_API}/entry/interpro/"
    f"{interpro_id}/"
)

entry_payload = request_json(
    session,
    entry_url
)

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
    family_name = (
        name_field
        or entry_meta.get("short_name")
        or interpro_id
    )

entry_type = (
    entry_meta.get("type")
    or entry_meta.get("entry_type")
)

print(
    f"InterPro název: {family_name}"
)

if entry_type:
    print(
        f"Typ entry: {entry_type}"
    )


# ------------------------------------------------------------
# TAXONOMICKÝ endpoint — žádné stahování všech proteinů
# ------------------------------------------------------------

initial_url = (
    f"{INTERPRO_API}/taxonomy/uniprot/"
    f"entry/interpro/{interpro_id}/"
)

url = initial_url

params = {
    "page_size": page_size,
    "extra_fields": (
        "scientific_name,full_name,rank,counters"
    ),
}

rows = []
seen_taxids = set()
page = 0
expected_count = None


def save_partial(rows_to_save):
    df = pd.DataFrame(rows_to_save)

    if df.empty:
        df = pd.DataFrame(
            columns=[
                "taxid",
                "organism",
                "rank",
                "family_member_count",
                "interpro_id",
                "interpro_name",
            ]
        )

    tmp = partial_file.with_suffix(
        partial_file.suffix + ".tmp"
    )

    df.to_csv(
        tmp,
        sep="\t",
        index=False
    )

    tmp.replace(partial_file)


def save_checkpoint(next_url, completed_page):
    payload = {
        "project_id": project_id,
        "interpro_id": interpro_id,
        "page_size": page_size,
        "completed_page": int(completed_page),
        "expected_taxa_count": expected_count,
        "next_url": next_url,
        "saved_taxa": len(rows),
    }

    tmp = checkpoint_file.with_suffix(
        checkpoint_file.suffix + ".tmp"
    )

    tmp.write_text(
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    tmp.replace(checkpoint_file)


# ------------------------------------------------------------
# Resume checkpoint
# ------------------------------------------------------------

if args.fresh:
    for stale in (
        checkpoint_file,
        partial_file,
    ):
        if stale.exists():
            stale.unlink()

elif (
    checkpoint_file.exists()
    and partial_file.exists()
):
    try:
        checkpoint = json.loads(
            checkpoint_file.read_text(
                encoding="utf-8"
            )
        )

        matches = (
            checkpoint.get("project_id")
            == project_id
            and checkpoint.get("interpro_id")
            == interpro_id
            and checkpoint.get("next_url")
        )

        if matches:
            partial_df = pd.read_csv(
                partial_file,
                sep="\t"
            )

            rows = (
                partial_df
                .where(
                    pd.notnull(partial_df),
                    None
                )
                .to_dict(
                    orient="records"
                )
            )

            seen_taxids = {
                int(row["taxid"])
                for row in rows
                if row.get("taxid") is not None
            }

            page = int(
                checkpoint.get(
                    "completed_page",
                    0
                )
            )

            expected_count = (
                checkpoint.get(
                    "expected_taxa_count"
                )
            )

            url = checkpoint.get(
                "next_url"
            )

            params = None

            print()
            print(
                "Nalezen rozpracovaný taxonomy download — "
                "pokračuji z checkpointu."
            )
            print(
                f"  hotových stránek: {page}"
            )
            print(
                f"  uložených taxonů: {len(rows)}"
            )
            print()

    except Exception as error:
        print(
            "Checkpoint se nepodařilo načíst; "
            f"stahuji od začátku ({error})."
        )

        rows = []
        seen_taxids = set()
        page = 0
        expected_count = None
        url = initial_url

        params = {
            "page_size": page_size,
            "extra_fields": (
                "scientific_name,"
                "full_name,"
                "rank,"
                "counters"
            ),
        }


# ------------------------------------------------------------
# Download taxonů
# ------------------------------------------------------------

while url:
    page += 1

    try:
        payload = request_json(
            session,
            url,
            params=params
        )

    except Exception:
        save_partial(rows)
        save_checkpoint(
            url,
            page - 1
        )

        print()
        print(
            "Taxonomy download byl přerušen, "
            "ale checkpoint je uložen."
        )
        print(
            f"Další spuštění naváže od stránky {page}."
        )
        print(
            f"Checkpoint: {checkpoint_file}"
        )
        raise

    # Od druhé stránky používáme next URL přesně tak,
    # jak ji vrátí InterPro API.
    params = None

    if expected_count is None:
        expected_count = payload.get("count")

        if expected_count is not None:
            print(
                f"InterPro hlásí "
                f"{expected_count} taxonomických záznamů."
            )

    results = (
        payload.get("results", [])
        or []
    )

    added = 0

    for item in results:
        taxid = parse_taxid(item)

        if (
            taxid is None
            or taxid in seen_taxids
        ):
            continue

        seen_taxids.add(taxid)
        added += 1

        protein_count = find_protein_count(
            item
        )

        rows.append({
            "taxid": int(taxid),
            "organism": (
                parse_scientific_name(item)
            ),
            "rank": parse_rank(item),
            "family_member_count": (
                int(protein_count)
                if protein_count is not None
                else None
            ),
            "interpro_id": interpro_id,
            "interpro_name": family_name,
        })

    next_url = payload.get("next")

    save_partial(rows)

    save_checkpoint(
        next_url,
        page
    )

    print(
        f"  stránka {page}: "
        f"+{added} taxonů; "
        f"unikátních taxonů {len(rows)} "
        "[checkpoint uložen]"
    )

    url = next_url

    if url:
        time.sleep(
            max(
                0.0,
                args.sleep
            )
        )


# ------------------------------------------------------------
# Finální TSV
# ------------------------------------------------------------

taxa_df = pd.DataFrame(rows)

if taxa_df.empty:
    taxa_df = pd.DataFrame(
        columns=[
            "taxid",
            "organism",
            "rank",
            "family_member_count",
            "interpro_id",
            "interpro_name",
        ]
    )

else:
    taxa_df["taxid"] = (
        taxa_df["taxid"]
        .astype(int)
    )

    taxa_df = taxa_df.sort_values(
        [
            "organism",
            "taxid",
        ],
        na_position="last"
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

    count_value = row.get(
        "family_member_count"
    )

    if pd.isna(count_value):
        count_value = None
    else:
        count_value = int(
            count_value
        )

    viewer_taxa[
        str(taxid)
    ] = {
        "taxid": taxid,
        "organism": (
            row.get("organism")
            or ""
        ),
        "rank": (
            row.get("rank")
            or ""
        ),
        "family_member": True,
        "family_member_count": count_value,
        "interpro_id": interpro_id,
        "interpro_name": family_name,

        # Taxon-first režim záměrně netahá jednotlivé UniProt ID.
        "uniprot_accessions": [],
    }


evidence = {
    "source": "InterPro taxonomy API",
    "mode": "taxonomy-first",
    "interpro_id": interpro_id,
    "interpro_name": family_name,
    "entry_type": entry_type,
    "taxa": viewer_taxa,
}

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
    "mode": "taxonomy-first",
    "api_reported_taxonomy_count": expected_count,
    "unique_taxids": int(
        taxa_df["taxid"].nunique()
    ) if not taxa_df.empty else 0,
    "files": {
        "taxa_tsv": (
            "family/family_taxa.tsv"
        ),
        "evidence_json": (
            "family/family_evidence.json"
        ),
    },
}

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
    "mode": "taxonomy-first",
    "interpro_id": interpro_id,
    "interpro_name": family_name,
    "entry_type": entry_type,
    "summary": (
        "family/family_summary.json"
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


# Úspěšný konec: checkpoint už nepotřebujeme.
for finished_file in (
    checkpoint_file,
    partial_file,
):
    if finished_file.exists():
        finished_file.unlink()


print()
print("=" * 70)
print("InterPro TAXON-FIRST evidence hotová")
print("=" * 70)
print(
    f"InterPro: {interpro_id} — {family_name}"
)
print(
    f"Taxonů:   {len(taxa_df)}"
)
print(
    f"Výstup:   {taxa_file}"
)
print(
    f"Evidence: {evidence_file}"
)
print("=" * 70)
