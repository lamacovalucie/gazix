#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re
import sys
import time
from io import StringIO

import pandas as pd
import requests


# ============================================================
# panther.py
#
# PANTHER enrichment vrstva pro Lucka-Structures.
#
# Workflow:
#   projectN
#     -> AlphaFind CSV
#     -> odhad query/self-hitu podle nejlepšího TM-score
#     -> UniProt FASTA
#     -> PANTHER TreeGrafter
#     -> PANTHER family
#     -> skutečný graft placement / subfamily
#     -> kontrola proti accessionu query, pokud už je v PANTHER stromu
#     -> family MSA + family tree
#     -> panther_members.tsv + NCBI Tax ID
#     -> panther_taxa.tsv
#     -> panther_evidence.json pro viewer
#
# Výstupy:
#   data/processed/projectN/panther/
#
# Poznámka:
#   Skript zatím NEMERGUJE PANTHER data přímo do viewer JSONu.
# ============================================================


# ------------------------------------------------------------
# ARGUMENTY
# ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=(
        "Zjistí PANTHER family/subfamily pro query protein projektu "
        "a stáhne family tree/MSA + tabulku členů rodiny."
    )
)

parser.add_argument(
    "project",
    nargs="?",
    help="Např. project3. Jinak se použije data/current_project.txt."
)

parser.add_argument(
    "--query",
    help="Ruční UniProt accession query proteinu, např. Q9C6B8."
)

parser.add_argument(
    "--query-taxid",
    type=int,
    help="Volitelný NCBI Tax ID query organismu."
)

parser.add_argument(
    "--min-self-tm",
    type=float,
    default=0.90,
    help=(
        "Minimální TM-score nejlepšího AlphaFind hitu, "
        "aby byl použit jako kandidát na self-hit. Default: 0.90."
    )
)

args = parser.parse_args()


# ------------------------------------------------------------
# CESTY
# ------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"

ALPHAFIND_ROOT_DIR = RAW_DIR / "alphafind"
PROCESSED_ROOT_DIR = DATA_DIR / "processed"
CURRENT_PROJECT_FILE = DATA_DIR / "current_project.txt"

PANTHER_API_ROOT = "https://www.pantherdb.org/services/oai/pantherdb"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

HEADERS = {
    "User-Agent": (
        "Lucka-Structures Panther enrichment "
        "(research/educational use)"
    ),
    "Accept": "application/json",
}


# ------------------------------------------------------------
# PROJECT
# ------------------------------------------------------------

PROJECT_RE = re.compile(r"^project(\d+)$", re.IGNORECASE)


def validate_project_id(project_id):
    project_id = str(project_id).strip().lower()

    if not PROJECT_RE.fullmatch(project_id):
        raise ValueError(
            f"Neplatné project ID: {project_id!r}. "
            "Použij project1, project2, ..."
        )

    return project_id


def get_project_id():
    if args.project:
        return validate_project_id(args.project)

    if not CURRENT_PROJECT_FILE.exists():
        raise FileNotFoundError(
            "Není zadán projectN a neexistuje data/current_project.txt."
        )

    value = CURRENT_PROJECT_FILE.read_text(
        encoding="utf-8"
    ).strip()

    return validate_project_id(value)


PROJECT_ID = get_project_id()

RAW_PROJECT_DIR = ALPHAFIND_ROOT_DIR / PROJECT_ID
PROCESSED_PROJECT_DIR = PROCESSED_ROOT_DIR / PROJECT_ID
PROJECT_JSON_FILE = PROCESSED_PROJECT_DIR / "project.json"
PANTHER_DIR = PROCESSED_PROJECT_DIR / "panther"

# Sdílená cache UniProt accession -> NCBI Tax ID.
# Je společná pro všechny projectN, aby se stejné accessiony
# nemusely dotazovat znovu při každém projektu.
UNIPROT_CACHE_DIR = RAW_DIR / "uniprot"
UNIPROT_TAXID_CACHE_FILE = (
    UNIPROT_CACHE_DIR / "uniprot_taxid_cache.tsv"
)

PANTHER_DIR.mkdir(parents=True, exist_ok=True)
UNIPROT_CACHE_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# PROJECT.JSON
# ------------------------------------------------------------

if PROJECT_JSON_FILE.exists():
    with open(PROJECT_JSON_FILE, "r", encoding="utf-8") as f:
        project_metadata = json.load(f)
else:
    project_metadata = {
        "project_id": PROJECT_ID
    }


# ------------------------------------------------------------
# ALPHAFIND CSV
# ------------------------------------------------------------

csv_files = sorted(RAW_PROJECT_DIR.glob("*.csv"))

if len(csv_files) == 0:
    raise FileNotFoundError(
        f"Ve {RAW_PROJECT_DIR} není žádný CSV."
    )

if len(csv_files) > 1:
    raise RuntimeError(
        f"Ve {RAW_PROJECT_DIR} je více CSV. "
        "Pro jeden projekt očekávám právě jeden."
    )

INPUT_CSV = csv_files[0]

df = pd.read_csv(INPUT_CSV, sep=";")


# ------------------------------------------------------------
# POMOCNÉ FUNKCE PRO SLOUPCE
# ------------------------------------------------------------

def normalize_column_name(value):
    return " ".join(
        str(value)
        .strip()
        .lower()
        .replace("_", " ")
        .replace("-", " ")
        .split()
    )


def find_column(candidates, required=True):
    normalized = {
        normalize_column_name(column): column
        for column in df.columns
    }

    for candidate in candidates:
        key = normalize_column_name(candidate)

        if key in normalized:
            return normalized[key]

    if required:
        raise KeyError(
            "Nenašla jsem očekávaný sloupec.\n"
            f"Hledala jsem: {candidates}\n"
            f"CSV obsahuje: {df.columns.tolist()}"
        )

    return None


TM_SCORE_COL = find_column([
    "TM-Score",
    "TM-score",
    "TM score",
    "tm_score",
])

UNIPROT_COL = find_column([
    "UniProt ID",
    "UniProt",
    "Protein ID",
    "Entry",
    "Accession",
])

TAXID_COL = find_column([
    "Tax ID",
    "TaxID",
    "Taxonomy ID",
    "Taxonomy",
], required=False)

ORGANISM_COL = find_column([
    "Organism",
    "Species",
], required=False)

RMSD_COL = find_column(
    ["RMSD"],
    required=False
)


def clean_numeric(series):
    return pd.to_numeric(
        series
        .astype(str)
        .str.replace("%", "", regex=False)
        .str.replace(",", ".", regex=False),
        errors="coerce"
    )


def normalize_uniprot(value):
    if pd.isna(value):
        return None

    text = str(value).strip()

    if not text:
        return None

    match = re.search(
        r"(?:sp|tr)\|([A-Z0-9]+)\|",
        text,
        re.IGNORECASE
    )

    if match:
        return match.group(1).upper()

    return text.split()[0].strip().upper()


df["_tm_numeric"] = clean_numeric(df[TM_SCORE_COL])


# ------------------------------------------------------------
# ODHAD QUERY / SELF-HITU
# ------------------------------------------------------------

def infer_query_from_alphafind():
    candidates = (
        df
        .dropna(subset=["_tm_numeric", UNIPROT_COL])
        .sort_values("_tm_numeric", ascending=False)
    )

    if candidates.empty:
        raise RuntimeError(
            "V CSV není použitelný TM-score + UniProt ID."
        )

    best = candidates.iloc[0]

    tm_score = float(best["_tm_numeric"])
    accession = normalize_uniprot(best[UNIPROT_COL])

    if accession is None or tm_score < args.min_self_tm:
        raise RuntimeError(
            "Automatická detekce query není dost jistá.\n"
            f"Nejlepší hit: {accession}, TM-score={tm_score:.3f}\n\n"
            "Zadej query ručně například:\n"
            "  python panther.py --query Q9C6B8"
        )

    taxid = None
    if TAXID_COL is not None:
        try:
            taxid = int(float(best[TAXID_COL]))
        except (TypeError, ValueError):
            pass

    organism = None
    if ORGANISM_COL is not None:
        value = best[ORGANISM_COL]
        if not pd.isna(value):
            organism = str(value)

    rmsd = None
    if RMSD_COL is not None:
        try:
            rmsd = float(str(best[RMSD_COL]).replace(",", "."))
        except (TypeError, ValueError):
            pass

    return {
        "uniprot": accession,
        "taxid": taxid,
        "organism": organism,
        "tm_score": tm_score,
        "rmsd": rmsd,
        "source": "alphafind_best_hit",
    }


# Priorita:
# 1) --query
# 2) query_uniprot v project.json
# 3) nejlepší AlphaFind hit

if args.query:
    query_info = {
        "uniprot": normalize_uniprot(args.query),
        "taxid": args.query_taxid,
        "organism": None,
        "tm_score": None,
        "rmsd": None,
        "source": "command_line",
    }

elif project_metadata.get("query_uniprot"):
    query_info = {
        "uniprot": normalize_uniprot(
            project_metadata["query_uniprot"]
        ),
        "taxid": project_metadata.get("query_taxid"),
        "organism": project_metadata.get("query_organism"),
        "tm_score": None,
        "rmsd": None,
        "source": "project_json",
    }

else:
    query_info = infer_query_from_alphafind()


if args.query_taxid is not None:
    query_info["taxid"] = args.query_taxid


print()
print("=" * 70)
print("PANTHER ENRICHMENT")
print("=" * 70)
print(f"Projekt:       {PROJECT_ID}")
print(f"AlphaFind CSV: {INPUT_CSV.name}")
print(f"Query UniProt: {query_info['uniprot']}")
print(f"Query Tax ID:  {query_info['taxid']}")
print(f"Organismus:    {query_info['organism']}")
print(f"Zdroj query:   {query_info['source']}")

if query_info["tm_score"] is not None:
    print(
        "Self-hit kandidát TM-score: "
        f"{query_info['tm_score']:.3f}"
    )

if query_info["rmsd"] is not None:
    print(
        "Self-hit kandidát RMSD: "
        f"{query_info['rmsd']:.3f}"
    )

print("=" * 70)


# ------------------------------------------------------------
# ULOŽIT QUERY DO PROJECT.JSON
# ------------------------------------------------------------

project_metadata["query_uniprot"] = query_info["uniprot"]

if query_info["taxid"] is not None:
    project_metadata["query_taxid"] = int(query_info["taxid"])

if query_info["organism"]:
    project_metadata["query_organism"] = query_info["organism"]

with open(PROJECT_JSON_FILE, "w", encoding="utf-8") as f:
    json.dump(
        project_metadata,
        f,
        ensure_ascii=False,
        indent=2
    )


# ------------------------------------------------------------
# UNIPROT FASTA
# ------------------------------------------------------------

accession = query_info["uniprot"]

fasta_url = UNIPROT_FASTA_URL.format(accession=accession)

print()
print(f"Stahuji sekvenci z UniProt: {accession}")

response = requests.get(
    fasta_url,
    timeout=30,
    headers={
        "User-Agent": HEADERS["User-Agent"]
    }
)

response.raise_for_status()

fasta_text = response.text.strip()

if not fasta_text.startswith(">"):
    raise RuntimeError(
        "UniProt nevrátil očekávaný FASTA formát."
    )

sequence = "".join(
    line.strip()
    for line in fasta_text.splitlines()
    if not line.startswith(">")
)

QUERY_FASTA_FILE = PANTHER_DIR / "query.fasta"

QUERY_FASTA_FILE.write_text(
    fasta_text + "\n",
    encoding="utf-8"
)

print(f"Délka query sekvence: {len(sequence)} aa")


# ------------------------------------------------------------
# PANTHER GENEINFO
# ------------------------------------------------------------

geneinfo_data = None

if query_info["taxid"] is not None:
    print()
    print(
        "Zkouším PANTHER geneinfo "
        "pro UniProt + Tax ID..."
    )

    try:
        geneinfo_response = requests.post(
            f"{PANTHER_API_ROOT}/geneinfo",
            data={
                "geneInputList": accession,
                "organism": str(query_info["taxid"]),
            },
            headers=HEADERS,
            timeout=60
        )

        geneinfo_response.raise_for_status()

        geneinfo_data = geneinfo_response.json()

        with open(
            PANTHER_DIR / "geneinfo.json",
            "w",
            encoding="utf-8"
        ) as f:
            json.dump(
                geneinfo_data,
                f,
                ensure_ascii=False,
                indent=2
            )

        print("PANTHER geneinfo uložen.")

    except Exception as error:
        print(
            "POZOR: geneinfo se nepodařilo získat."
        )
        print(f"  {error}")
        print(
            "Pokračuji přes TreeGrafter."
        )


# ------------------------------------------------------------
# PANTHER TREEGRAFTER
# ------------------------------------------------------------

print()
print(
    "Posílám query sekvenci do "
    "PANTHER TreeGrafter..."
)

graft_response = requests.post(
    f"{PANTHER_API_ROOT}/graftsequence",
    data={
        "sequence": sequence
    },
    headers=HEADERS,
    timeout=180
)

graft_response.raise_for_status()

graft_data = graft_response.json()

GRAFT_FILE = PANTHER_DIR / "graftsequence.json"

with open(GRAFT_FILE, "w", encoding="utf-8") as f:
    json.dump(
        graft_data,
        f,
        ensure_ascii=False,
        indent=2
    )

print("TreeGrafter odpověď uložena.")


# ------------------------------------------------------------
# POMOCNÉ FUNKCE PRO PANTHER TREE
# ------------------------------------------------------------

PANTHER_FAMILY_RE = re.compile(
    r"\b(PTHR\d+)\b",
    re.IGNORECASE
)

UNIPROT_IN_NODE_RE = re.compile(
    r"(?:^|\|)UniProtKB=([A-Z0-9]+)(?:\||$)",
    re.IGNORECASE
)


def children_of(node):
    """
    Vrátí seznam potomků PANTHER annotation_node.
    PANTHER JSON používá:
      node["children"]["annotation_node"]
    ale někdy může být jediný potomek reprezentován dict místo list.
    """
    if not isinstance(node, dict):
        return []

    children = node.get("children")

    if not isinstance(children, dict):
        return []

    annotation = children.get("annotation_node", [])

    if isinstance(annotation, list):
        return annotation

    if isinstance(annotation, dict):
        return [annotation]

    return []


def walk_nodes(obj):
    """
    Rekurzivně projde všechny dicty v JSONu.
    """
    if isinstance(obj, dict):
        yield obj

        for value in obj.values():
            yield from walk_nodes(value)

    elif isinstance(obj, list):
        for value in obj:
            yield from walk_nodes(value)


def is_grafted_node(node):
    if not isinstance(node, dict):
        return False

    node_name = str(node.get("node_name", ""))
    accession_value = str(node.get("accession", ""))
    public_id = str(node.get("public_id", ""))

    return (
        node_name.startswith("GRAFTED|")
        or accession_value == "ANGRAFTED"
        or public_id == "PTNGRAFTED"
    )


def extract_uniprot_from_node(node):
    if not isinstance(node, dict):
        return None

    node_name = str(node.get("node_name", ""))

    match = UNIPROT_IN_NODE_RE.search(node_name)

    if match:
        return match.group(1).upper()

    return None


def find_query_accession_node(tree_data, query_accession):
    """
    Najde existující PANTHER leaf se stejným UniProt accessionem,
    pokud query protein už v rodinném stromu je.
    """
    query_accession = query_accession.upper()

    for node in walk_nodes(tree_data):
        node_uniprot = extract_uniprot_from_node(node)

        if node_uniprot == query_accession:
            return node

    return None


def collect_subfamily_ids_from_subtree(node):
    """
    Sesbírá sf_id / prop_sf_id z daného podstromu.
    Používá se jen lokálně kolem graft pozice,
    nikoli přes celý family tree.
    """
    ids = set()

    for item in walk_nodes(node):
        for key in ("sf_id", "prop_sf_id"):
            value = item.get(key)

            if (
                isinstance(value, str)
                and value.upper().startswith("PTHR")
                and ":SF" in value.upper()
            ):
                ids.add(value.upper())

    return ids


def find_graft_placement(tree_data):
    """
    Najde rodičovský uzel, pod kterým je vložen GRAFTED node.

    Primární strategie:
      1) pokud rodič sám nese sf_id/prop_sf_id, použij ho
      2) jinak prohlédni NE-GRAFTED sourozence
         a hledej nejbližší subfamily ID
      3) pokud všichni sourozenci ukazují na jedinou subfamily,
         ber ji jako graft placement

    Pro PIN1 ukazuje sibling větev na PTHR31752:SF18.
    """
    placement = {
        "subfamily_id": None,
        "subfamily_name": None,
        "method": None,
        "parent_node": None,
        "sibling_candidates": [],
    }

    def recurse(node):
        if not isinstance(node, dict):
            return None

        children = children_of(node)

        if children:
            grafted_children = [
                child
                for child in children
                if is_grafted_node(child)
            ]

            if grafted_children:
                # 1) rodič má vlastní subfamily anotaci
                parent_sf = (
                    node.get("sf_id")
                    or node.get("prop_sf_id")
                )

                if isinstance(parent_sf, str) and ":SF" in parent_sf.upper():
                    placement["subfamily_id"] = parent_sf.upper()
                    placement["subfamily_name"] = node.get("sf_name")
                    placement["method"] = "grafted_parent"
                    placement["parent_node"] = node
                    return placement

                # 2) podíváme se na bezprostřední sourozence
                siblings = [
                    child
                    for child in children
                    if not is_grafted_node(child)
                ]

                direct_candidates = []

                for sibling in siblings:
                    sf = (
                        sibling.get("sf_id")
                        or sibling.get("prop_sf_id")
                    )

                    if isinstance(sf, str) and ":SF" in sf.upper():
                        direct_candidates.append({
                            "subfamily_id": sf.upper(),
                            "subfamily_name": sibling.get("sf_name"),
                            "source": "direct_sibling",
                        })
                    else:
                        subtree_ids = sorted(
                            collect_subfamily_ids_from_subtree(sibling)
                        )

                        for subtree_id in subtree_ids:
                            direct_candidates.append({
                                "subfamily_id": subtree_id,
                                "subfamily_name": sibling.get("sf_name"),
                                "source": "sibling_subtree",
                            })

                placement["sibling_candidates"] = direct_candidates

                unique_ids = sorted({
                    item["subfamily_id"]
                    for item in direct_candidates
                })

                if len(unique_ids) == 1:
                    selected = unique_ids[0]

                    names = [
                        item["subfamily_name"]
                        for item in direct_candidates
                        if (
                            item["subfamily_id"] == selected
                            and item["subfamily_name"]
                        )
                    ]

                    placement["subfamily_id"] = selected
                    placement["subfamily_name"] = (
                        names[0]
                        if names
                        else None
                    )
                    placement["method"] = "grafted_sibling"
                    placement["parent_node"] = node

                else:
                    placement["method"] = "grafted_ambiguous"
                    placement["parent_node"] = node

                return placement

            for child in children:
                found = recurse(child)

                if found is not None:
                    return found

        return None

    # TreeGrafter odpověď má typicky hlavní strom pod search.annotation_node
    search = tree_data.get("search", {}) if isinstance(tree_data, dict) else {}
    root = search.get("annotation_node")

    if root is None:
        root = tree_data

    result = recurse(root)

    if result is None:
        placement["method"] = "grafted_not_found"

    return placement


def family_id_from_graft(tree_data):
    """
    Preferuje search.book, např. PTHR31752.
    """
    if isinstance(tree_data, dict):
        search = tree_data.get("search", {})

        if isinstance(search, dict):
            book = search.get("book")

            if isinstance(book, str):
                match = PANTHER_FAMILY_RE.search(book)

                if match:
                    return match.group(1).upper()

    # fallback: první PTHR family ID v JSONu
    for node in walk_nodes(tree_data):
        for value in node.values():
            if isinstance(value, str):
                match = PANTHER_FAMILY_RE.search(value)

                if match:
                    return match.group(1).upper()

    return None


def family_name_from_graft(tree_data):
    if not isinstance(tree_data, dict):
        return None

    search = tree_data.get("search", {})

    if not isinstance(search, dict):
        return None

    return search.get("family_name")


# ------------------------------------------------------------
# FAMILY + GRAFT PLACEMENT
# ------------------------------------------------------------

PANTHER_FAMILY_ID = family_id_from_graft(graft_data)

if not PANTHER_FAMILY_ID:
    print()
    print(
        "Nepodařilo se automaticky najít PANTHER family ID."
    )
    print("Raw odpověď je uložená zde:")
    print(f"  {GRAFT_FILE}")
    sys.exit(0)


PANTHER_FAMILY_NAME = family_name_from_graft(graft_data)

graft_placement = find_graft_placement(graft_data)

GRAFT_SUBFAMILY_ID = graft_placement["subfamily_id"]
GRAFT_SUBFAMILY_NAME = graft_placement["subfamily_name"]


# ------------------------------------------------------------
# VALIDACE PODLE ACCESSIONU QUERY
# ------------------------------------------------------------

query_tree_node = find_query_accession_node(
    graft_data,
    accession
)

ACCESSION_SUBFAMILY_ID = None
ACCESSION_SUBFAMILY_NAME = None
ACCESSION_VALIDATION = "query_not_present_in_panther_tree"

if query_tree_node is not None:
    accession_sf = (
        query_tree_node.get("prop_sf_id")
        or query_tree_node.get("sf_id")
    )

    if isinstance(accession_sf, str) and ":SF" in accession_sf.upper():
        ACCESSION_SUBFAMILY_ID = accession_sf.upper()

    ACCESSION_SUBFAMILY_NAME = query_tree_node.get("sf_name")

    if (
        GRAFT_SUBFAMILY_ID
        and ACCESSION_SUBFAMILY_ID
        and GRAFT_SUBFAMILY_ID == ACCESSION_SUBFAMILY_ID
    ):
        ACCESSION_VALIDATION = "graft_and_accession_agree"

    elif (
        GRAFT_SUBFAMILY_ID
        and ACCESSION_SUBFAMILY_ID
        and GRAFT_SUBFAMILY_ID != ACCESSION_SUBFAMILY_ID
    ):
        ACCESSION_VALIDATION = "graft_and_accession_disagree"

    else:
        ACCESSION_VALIDATION = "accession_found_but_graft_unresolved"


# Pokud graft placement selže, ale accession query už v PANTHER stromu je,
# použijeme accession jako fallback. V summary ale zachováme, že šlo o fallback.

FINAL_SUBFAMILY_ID = GRAFT_SUBFAMILY_ID
FINAL_SUBFAMILY_NAME = GRAFT_SUBFAMILY_NAME
SUBFAMILY_SOURCE = graft_placement["method"]

if FINAL_SUBFAMILY_ID is None and ACCESSION_SUBFAMILY_ID is not None:
    FINAL_SUBFAMILY_ID = ACCESSION_SUBFAMILY_ID
    FINAL_SUBFAMILY_NAME = ACCESSION_SUBFAMILY_NAME
    SUBFAMILY_SOURCE = "query_accession_fallback"

if (
    FINAL_SUBFAMILY_NAME is None
    and ACCESSION_SUBFAMILY_ID == FINAL_SUBFAMILY_ID
):
    FINAL_SUBFAMILY_NAME = ACCESSION_SUBFAMILY_NAME


print()
print("PANTHER zařazení query:")
print(f"  Family:             {PANTHER_FAMILY_ID}")
print(f"  Family name:        {PANTHER_FAMILY_NAME}")
print(f"  Graft subfamily:    {GRAFT_SUBFAMILY_ID}")
print(f"  Graft method:       {graft_placement['method']}")
print(f"  Accession subfamily:{ACCESSION_SUBFAMILY_ID}")
print(f"  Kontrola:           {ACCESSION_VALIDATION}")
print(f"  Výsledná subfamily: {FINAL_SUBFAMILY_ID}")
print(f"  Subfamily name:     {FINAL_SUBFAMILY_NAME}")

if ACCESSION_VALIDATION == "graft_and_accession_disagree":
    print()
    print(
        "POZOR: graft placement a existující accession v PANTHER stromu "
        "ukazují na různé subfamilies."
    )
    print(
        "Výsledek automaticky nepřepisuji podle accessionu; "
        "zkontroluj graftsequence.json."
    )


# ------------------------------------------------------------
# FAMILY MSA + TREE
# ------------------------------------------------------------

def panther_get(endpoint, params):
    response = requests.get(
        f"{PANTHER_API_ROOT}/{endpoint}",
        params=params,
        headers=HEADERS,
        timeout=180
    )

    response.raise_for_status()
    return response.json()


print()
print(
    "Stahuji PANTHER family MSA: "
    f"{PANTHER_FAMILY_ID}"
)

msa_data = panther_get(
    "familymsa",
    {
        "family": PANTHER_FAMILY_ID
    }
)

with open(
    PANTHER_DIR / "family_msa.json",
    "w",
    encoding="utf-8"
) as f:
    json.dump(
        msa_data,
        f,
        ensure_ascii=False,
        indent=2
    )


print(
    "Stahuji PANTHER family tree: "
    f"{PANTHER_FAMILY_ID}"
)

tree_data = panther_get(
    "treeinfo",
    {
        "family": PANTHER_FAMILY_ID
    }
)

with open(
    PANTHER_DIR / "family_tree.json",
    "w",
    encoding="utf-8"
) as f:
    json.dump(
        tree_data,
        f,
        ensure_ascii=False,
        indent=2
    )


# ------------------------------------------------------------
# PANTHER MEMBERS TSV
# ------------------------------------------------------------

def leaf_nodes(tree_data):
    """
    Vrátí biologické leaf nodes; vynechá syntetický GRAFTED node.
    """
    leaves = []

    for node in walk_nodes(tree_data):
        if is_grafted_node(node):
            continue

        if node.get("tree_node_type") == "LEAF":
            leaves.append(node)

    return leaves


def make_members_dataframe(tree_data, query_accession):
    rows = []

    for node in leaf_nodes(tree_data):
        uniprot = extract_uniprot_from_node(node)

        subfamily_id = (
            node.get("prop_sf_id")
            or node.get("sf_id")
        )

        if isinstance(subfamily_id, str):
            subfamily_id = subfamily_id.upper()

        rows.append({
            "uniprot": uniprot,
            "organism": node.get("organism"),
            "gene_symbol": node.get("gene_symbol"),
            "definition": node.get("definition"),
            "subfamily_id": subfamily_id,
            "subfamily_name": node.get("sf_name"),
            "gene_id": node.get("gene_id"),
            "panther_accession": node.get("accession"),
            "panther_public_id": node.get("public_id"),
            "branch_length": node.get("branch_length"),
            "is_query_accession": (
                bool(uniprot)
                and uniprot.upper() == query_accession.upper()
            ),
        })

    members_df = pd.DataFrame(rows)

    if not members_df.empty:
        members_df = members_df.sort_values(
            by=[
                "subfamily_id",
                "organism",
                "uniprot",
            ],
            na_position="last"
        ).reset_index(drop=True)

    return members_df


# Pro členy používáme graft_data, protože v našem ověřeném PANTHER
# TreeGrafter JSONu už je celý rodinný strom včetně sf_id/prop_sf_id,
# organismů, gene symbols a UniProt accessionů.
members_df = make_members_dataframe(
    graft_data,
    accession
)

MEMBERS_FILE = PANTHER_DIR / "panther_members.tsv"

members_df.to_csv(
    MEMBERS_FILE,
    sep="\t",
    index=False
)

unique_organisms = (
    members_df["organism"].dropna().nunique()
    if "organism" in members_df.columns
    else 0
)

unique_subfamilies = (
    members_df["subfamily_id"].dropna().nunique()
    if "subfamily_id" in members_df.columns
    else 0
)

query_rows_count = (
    int(members_df["is_query_accession"].sum())
    if "is_query_accession" in members_df.columns
    else 0
)

print()
print("PANTHER members TSV:")
print(f"  Členů rodiny:       {len(members_df)}")
print(f"  Organismů:          {unique_organisms}")
print(f"  Subfamilies:        {unique_subfamilies}")
print(f"  Query accession řádků: {query_rows_count}")
print(f"  Soubor:             {MEMBERS_FILE}")


# ------------------------------------------------------------
# DOPLNIT NCBI TAX ID PŘES UNIPROT
# ------------------------------------------------------------

def load_uniprot_taxid_cache():
    """
    Načte sdílenou cache:
        uniprot  taxid  organism_uniprot
    """
    if not UNIPROT_TAXID_CACHE_FILE.exists():
        return pd.DataFrame(
            columns=[
                "uniprot",
                "taxid",
                "organism_uniprot",
            ]
        )

    try:
        cache = pd.read_csv(
            UNIPROT_TAXID_CACHE_FILE,
            sep="\t",
            dtype={
                "uniprot": "string",
                "taxid": "Int64",
                "organism_uniprot": "string",
            }
        )
    except Exception as error:
        print()
        print(
            "POZOR: Nepodařilo se načíst UniProt taxID cache."
        )
        print(f"  {error}")
        print("Pokračuji s prázdnou cache.")
        return pd.DataFrame(
            columns=[
                "uniprot",
                "taxid",
                "organism_uniprot",
            ]
        )

    for column in [
        "uniprot",
        "taxid",
        "organism_uniprot",
    ]:
        if column not in cache.columns:
            cache[column] = pd.NA

    cache["uniprot"] = (
        cache["uniprot"]
        .astype("string")
        .str.upper()
    )

    cache["taxid"] = pd.to_numeric(
        cache["taxid"],
        errors="coerce"
    ).astype("Int64")

    return cache[
        [
            "uniprot",
            "taxid",
            "organism_uniprot",
        ]
    ]


def save_uniprot_taxid_cache(cache):
    """
    Atomický zápis cache, aby se při přerušení nepoškodila.
    """
    cache = cache.copy()

    cache["uniprot"] = (
        cache["uniprot"]
        .astype("string")
        .str.upper()
    )

    cache["taxid"] = pd.to_numeric(
        cache["taxid"],
        errors="coerce"
    ).astype("Int64")

    cache = (
        cache
        .dropna(subset=["uniprot"])
        .drop_duplicates(
            subset=["uniprot"],
            keep="last"
        )
        .sort_values("uniprot")
        .reset_index(drop=True)
    )

    temp_file = (
        UNIPROT_TAXID_CACHE_FILE
        .with_suffix(".tmp.tsv")
    )

    cache.to_csv(
        temp_file,
        sep="\t",
        index=False
    )

    temp_file.replace(
        UNIPROT_TAXID_CACHE_FILE
    )


def chunks(values, size):
    for start in range(
        0,
        len(values),
        size
    ):
        yield values[
            start:start + size
        ]


def fetch_uniprot_taxids(accessions):
    """
    Dávkově stáhne accession -> organism_id z UniProt REST API.

    Používá search endpoint místo 1 request / 1 protein.
    """
    accessions = sorted({
        str(value).strip().upper()
        for value in accessions
        if (
            value is not None
            and not pd.isna(value)
            and str(value).strip()
        )
    })

    if not accessions:
        return pd.DataFrame(
            columns=[
                "uniprot",
                "taxid",
                "organism_uniprot",
            ]
        )

    result_frames = []

    # Menší dávky = kratší URL a robustnější requesty.
    batch_size = 50

    for batch_number, batch in enumerate(
        chunks(
            accessions,
            batch_size
        ),
        start=1
    ):
        query = " OR ".join(
            f"accession:{acc}"
            for acc in batch
        )

        params = {
            "query": f"({query})",
            "format": "tsv",
            "fields": (
                "accession,"
                "organism_id,"
                "organism_name"
            ),
            "size": len(batch),
        }

        url = (
            "https://rest.uniprot.org/"
            "uniprotkb/search"
        )

        print(
            f"  UniProt taxID batch "
            f"{batch_number}: {len(batch)} accessionů"
        )

        last_error = None

        for attempt in range(1, 4):
            try:
                response = requests.get(
                    url,
                    params=params,
                    headers={
                        "User-Agent":
                            HEADERS["User-Agent"]
                    },
                    timeout=90
                )

                if response.status_code == 429:
                    retry_after = (
                        response.headers
                        .get("Retry-After")
                    )

                    try:
                        wait_seconds = min(
                            int(retry_after),
                            60
                        )
                    except (
                        TypeError,
                        ValueError
                    ):
                        wait_seconds = 10 * attempt

                    print(
                        "    UniProt 429; "
                        f"čekám {wait_seconds} s..."
                    )

                    time.sleep(
                        wait_seconds
                    )
                    continue

                response.raise_for_status()

                batch_df = pd.read_csv(
                    StringIO(
                        response.text
                    ),
                    sep="\t"
                )

                # UniProt TSV názvy sloupců:
                # Entry, Organism (ID), Organism
                rename_map = {
                    "Entry": "uniprot",
                    "Organism (ID)": "taxid",
                    "Organism": "organism_uniprot",
                }

                batch_df = batch_df.rename(
                    columns=rename_map
                )

                for column in [
                    "uniprot",
                    "taxid",
                    "organism_uniprot",
                ]:
                    if column not in batch_df.columns:
                        batch_df[column] = pd.NA

                batch_df["uniprot"] = (
                    batch_df["uniprot"]
                    .astype("string")
                    .str.upper()
                )

                batch_df["taxid"] = (
                    pd.to_numeric(
                        batch_df["taxid"],
                        errors="coerce"
                    )
                    .astype("Int64")
                )

                result_frames.append(
                    batch_df[
                        [
                            "uniprot",
                            "taxid",
                            "organism_uniprot",
                        ]
                    ]
                )

                last_error = None
                break

            except Exception as error:
                last_error = error

                if attempt < 3:
                    wait_seconds = (
                        5 * attempt
                    )

                    print(
                        "    request selhal; "
                        f"další pokus za "
                        f"{wait_seconds} s"
                    )

                    time.sleep(
                        wait_seconds
                    )

        if last_error is not None:
            print(
                "    POZOR: tento batch "
                "se nepodařilo stáhnout."
            )
            print(
                f"    {last_error}"
            )

        # Krátká pauza mezi dávkami.
        time.sleep(0.3)

    if not result_frames:
        return pd.DataFrame(
            columns=[
                "uniprot",
                "taxid",
                "organism_uniprot",
            ]
        )

    result = pd.concat(
        result_frames,
        ignore_index=True
    )

    result = (
        result
        .dropna(subset=["uniprot"])
        .drop_duplicates(
            subset=["uniprot"],
            keep="last"
        )
    )

    return result


print()
print("Doplňuji NCBI Tax ID členům PANTHER rodiny...")

uniprot_cache = (
    load_uniprot_taxid_cache()
)

cached_map = {
    str(row["uniprot"]).upper():
        (
            row["taxid"],
            row["organism_uniprot"],
        )
    for _, row in uniprot_cache.iterrows()
    if not pd.isna(
        row["uniprot"]
    )
}

member_accessions = sorted({
    str(value).upper()
    for value in members_df["uniprot"].dropna()
    if str(value).strip()
})

missing_accessions = [
    acc
    for acc in member_accessions
    if acc not in cached_map
]

print(
    f"  UniProt accessionů v TSV: "
    f"{len(member_accessions)}"
)
print(
    f"  Už v cache: "
    f"{len(member_accessions) - len(missing_accessions)}"
)
print(
    f"  Nově k dotazu: "
    f"{len(missing_accessions)}"
)

if missing_accessions:
    new_taxids = fetch_uniprot_taxids(
        missing_accessions
    )

    if not new_taxids.empty:
        uniprot_cache = pd.concat(
            [
                uniprot_cache,
                new_taxids,
            ],
            ignore_index=True
        )

        save_uniprot_taxid_cache(
            uniprot_cache
        )

        print(
            "  UniProt taxID cache "
            "aktualizována."
        )

# Znovu načíst finální cache po případném zápisu.
uniprot_cache = (
    load_uniprot_taxid_cache()
)

members_df = members_df.merge(
    uniprot_cache,
    on="uniprot",
    how="left"
)

members_df["taxid"] = pd.to_numeric(
    members_df["taxid"],
    errors="coerce"
).astype("Int64")

members_df[
    "same_query_subfamily"
] = (
    members_df["subfamily_id"]
    .astype("string")
    == str(
        FINAL_SUBFAMILY_ID
    )
)

# Přepsat panther_members.tsv už obohacenou verzí.
members_df.to_csv(
    MEMBERS_FILE,
    sep="\t",
    index=False
)

members_with_taxid = int(
    members_df["taxid"]
    .notna()
    .sum()
)

members_without_taxid = (
    len(members_df)
    - members_with_taxid
)

unique_taxids = int(
    members_df["taxid"]
    .dropna()
    .nunique()
)

print(
    f"  Tax ID doplněno:   "
    f"{members_with_taxid}/{len(members_df)}"
)
print(
    f"  Bez Tax ID:        "
    f"{members_without_taxid}"
)
print(
    f"  Unikátních Tax ID: "
    f"{unique_taxids}"
)


# ------------------------------------------------------------
# AGREGACE PRO VIEWER: panther_taxa.tsv
# ------------------------------------------------------------

def join_unique(values):
    cleaned = sorted({
        str(value)
        for value in values
        if (
            value is not None
            and not pd.isna(value)
            and str(value).strip()
        )
    })

    return ";".join(cleaned)


taxon_rows = []

members_with_known_taxid_df = (
    members_df[
        members_df["taxid"].notna()
    ]
    .copy()
)

for taxid, group in (
    members_with_known_taxid_df
    .groupby("taxid")
):
    organisms = [
        value
        for value in group[
            "organism"
        ].dropna()
    ]

    organism_name = (
        str(organisms[0])
        if organisms
        else None
    )

    family_count = len(group)

    query_sf_group = group[
        group[
            "same_query_subfamily"
        ]
    ]

    query_subfamily_count = (
        len(query_sf_group)
    )

    taxon_rows.append({
        "taxid": int(taxid),
        "organism": organism_name,
        "panther_family_id":
            PANTHER_FAMILY_ID,
        "family_member_count":
            family_count,
        "query_subfamily_id":
            FINAL_SUBFAMILY_ID,
        "query_subfamily_member_count":
            query_subfamily_count,
        "same_query_subfamily":
            query_subfamily_count > 0,
        "subfamily_ids":
            join_unique(
                group[
                    "subfamily_id"
                ]
            ),
        "uniprot_accessions":
            join_unique(
                group[
                    "uniprot"
                ]
            ),
    })


panther_taxa_df = pd.DataFrame(
    taxon_rows
)

if not panther_taxa_df.empty:
    panther_taxa_df = (
        panther_taxa_df
        .sort_values(
            [
                "same_query_subfamily",
                "family_member_count",
                "organism",
            ],
            ascending=[
                False,
                False,
                True,
            ]
        )
        .reset_index(drop=True)
    )

PANTHER_TAXA_FILE = (
    PANTHER_DIR
    / "panther_taxa.tsv"
)

panther_taxa_df.to_csv(
    PANTHER_TAXA_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# EVIDENCE JSON PRO VIEWER
# ------------------------------------------------------------

taxa_json = {}

for _, row in (
    panther_taxa_df.iterrows()
):
    taxid_key = str(
        int(row["taxid"])
    )

    taxa_json[
        taxid_key
    ] = {
        "taxid":
            int(row["taxid"]),

        "organism":
            (
                None
                if pd.isna(
                    row["organism"]
                )
                else str(
                    row["organism"]
                )
            ),

        "family_member":
            True,

        "family_member_count":
            int(
                row[
                    "family_member_count"
                ]
            ),

        "same_query_subfamily":
            bool(
                row[
                    "same_query_subfamily"
                ]
            ),

        "query_subfamily_member_count":
            int(
                row[
                    "query_subfamily_member_count"
                ]
            ),

        "subfamily_ids": [
            value
            for value in str(
                row[
                    "subfamily_ids"
                ]
            ).split(";")
            if value
        ],

        "uniprot_accessions": [
            value
            for value in str(
                row[
                    "uniprot_accessions"
                ]
            ).split(";")
            if value
        ],
    }


panther_evidence = {
    "source":
        "PANTHER",

    "project_id":
        PROJECT_ID,

    "query_uniprot":
        accession,

    "family": {
        "id":
            PANTHER_FAMILY_ID,
        "name":
            PANTHER_FAMILY_NAME,
    },

    "query_subfamily": {
        "id":
            FINAL_SUBFAMILY_ID,
        "name":
            FINAL_SUBFAMILY_NAME,
    },

    "legend_semantics": {
        "family_member":
            "light_blue",
        "same_query_subfamily":
            "dark_blue",
        "alphafind_and_panther":
            "purple",
    },

    "coverage_note": (
        "PANTHER gene-tree species coverage is "
        "reference-set based and is not an exhaustive "
        "survey of all species containing homologs."
    ),

    "taxa":
        taxa_json,
}


PANTHER_EVIDENCE_FILE = (
    PANTHER_DIR
    / "panther_evidence.json"
)

with open(
    PANTHER_EVIDENCE_FILE,
    "w",
    encoding="utf-8"
) as f:
    json.dump(
        panther_evidence,
        f,
        ensure_ascii=False,
        indent=2
    )


print()
print("PANTHER taxonomická evidence:")
print(
    f"  Taxonů pro viewer:  "
    f"{len(panther_taxa_df)}"
)
print(
    f"  TSV:                "
    f"{PANTHER_TAXA_FILE}"
)
print(
    f"  JSON:               "
    f"{PANTHER_EVIDENCE_FILE}"
)


# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

summary = {
    "project_id": PROJECT_ID,
    "query_uniprot": accession,
    "query_taxid": query_info["taxid"],
    "query_organism": query_info["organism"],
    "query_source": query_info["source"],
    "query_sequence_length": len(sequence),

    "panther_family_id": PANTHER_FAMILY_ID,
    "panther_family_name": PANTHER_FAMILY_NAME,

    "panther_subfamily_id": FINAL_SUBFAMILY_ID,
    "panther_subfamily_name": FINAL_SUBFAMILY_NAME,
    "panther_subfamily_source": SUBFAMILY_SOURCE,

    "graft": {
        "subfamily_id": GRAFT_SUBFAMILY_ID,
        "subfamily_name": GRAFT_SUBFAMILY_NAME,
        "method": graft_placement["method"],
        "sibling_candidates": graft_placement["sibling_candidates"],
    },

    "query_accession_validation": {
        "query_found_in_panther_tree": query_tree_node is not None,
        "accession_subfamily_id": ACCESSION_SUBFAMILY_ID,
        "accession_subfamily_name": ACCESSION_SUBFAMILY_NAME,
        "status": ACCESSION_VALIDATION,
    },

    "members": {
        "count": len(members_df),
        "unique_organisms": int(unique_organisms),
        "unique_subfamilies": int(unique_subfamilies),
        "query_accession_rows": int(query_rows_count),
        "with_taxid": int(members_with_taxid),
        "without_taxid": int(members_without_taxid),
        "unique_taxids": int(unique_taxids),
        "viewer_taxa": int(len(panther_taxa_df)),
    },

    "files": {
        "query_fasta": "query.fasta",
        "geneinfo": (
            "geneinfo.json"
            if geneinfo_data is not None
            else None
        ),
        "graftsequence": "graftsequence.json",
        "family_msa": "family_msa.json",
        "family_tree": "family_tree.json",
        "members_tsv": "panther_members.tsv",
        "taxa_tsv": "panther_taxa.tsv",
        "evidence_json": "panther_evidence.json",
    }
}

SUMMARY_FILE = PANTHER_DIR / "panther_summary.json"

with open(SUMMARY_FILE, "w", encoding="utf-8") as f:
    json.dump(
        summary,
        f,
        ensure_ascii=False,
        indent=2
    )


# ------------------------------------------------------------
# PROJECT.JSON
# ------------------------------------------------------------

project_metadata["panther"] = {
    "family_id": PANTHER_FAMILY_ID,
    "family_name": PANTHER_FAMILY_NAME,
    "subfamily_id": FINAL_SUBFAMILY_ID,
    "subfamily_name": FINAL_SUBFAMILY_NAME,
    "subfamily_source": SUBFAMILY_SOURCE,
    "validation": ACCESSION_VALIDATION,
    "summary": "panther/panther_summary.json",
    "members_tsv": "panther/panther_members.tsv",
    "taxa_tsv": "panther/panther_taxa.tsv",
    "evidence_json": "panther/panther_evidence.json",
}

with open(PROJECT_JSON_FILE, "w", encoding="utf-8") as f:
    json.dump(
        project_metadata,
        f,
        ensure_ascii=False,
        indent=2
    )


print()
print("=" * 70)
print("HOTOVO")
print("=" * 70)
print(f"Query UniProt:       {accession}")
print(f"PANTHER family:      {PANTHER_FAMILY_ID}")
print(f"Family name:         {PANTHER_FAMILY_NAME}")
print(f"PANTHER subfamily:   {FINAL_SUBFAMILY_ID}")
print(f"Subfamily name:      {FINAL_SUBFAMILY_NAME}")
print(f"Subfamily source:    {SUBFAMILY_SOURCE}")
print(f"Validace accessionu: {ACCESSION_VALIDATION}")
print(f"Členů v TSV:         {len(members_df)}")
print(f"Výstupy:             {PANTHER_DIR}")
print()
print(
    "Další krok: napojit panther_evidence.json na existující "
    "NCBI strom ve vieweru a kombinovat PANTHER + AlphaFind barvy."
)
