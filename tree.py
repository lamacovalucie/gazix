# skript z 15/9, resi vylepseni ukladani dat tak, aby kazda rodina proteinu mela svou slozku
#data/
#├── raw/
#│   ├── alphafind/
#│   │   ├── betv1/
#│   │   │   └── něco.csv
#│   │   └── 1fas/
#│   │       └── něco.csv
#│   ├── ncbi_taxonomy/       ← sdílené
#│   └── taxon_names_cs/      ← sdílené
#│
#└── processed/
#    ├── betv1/
#    │   └── všechny výstupy
#    └── 1fas/
#        └── všechny výstupy




# ==== import hotových knhoven, funkcí, balíčků atd. ===
from pathlib import Path #path pro práci s chytrými cestami; json je textový formát - kvůli ukládání metadat z OpenTree)
import json #knihovna pro práci se json
from datetime import datetime #z knihovny datetime si načti objekt datetime - převádí na datum a čas (kdy byl strom stažen)
import argparse #pro předání názvu projektu / režimu spuštění
import shutil #pro přesunutí nového AlphaFind CSV z inboxu do projectN
# pathlib,json a datetime jsou součástí balíčku pythonu, requests a pandas jsou externí co musí být instalovány v mém prosředí.

import requests # = načti knihovnu requests, která umí posílat požadavky přes internet (komunikace s OPENTREE přes API)
import pandas as pd # = načti knihovnu pandas a dej jí zkratku pd.
import tarfile #kvůli stažení NCBI taxonomie
import re
import math
import time


# ============================================================
# 1. CESTY V PROJEKTU + SPRÁVA PROJEKTŮ
# ============================================================

# PRINCIP:
#
# - technické ID projektu je stabilní:
#       project1, project2, project3, ...
#
# - lidský název projektu NENÍ název složky.
#   Je uložený v:
#       data/processed/projectN/project.json
#
# - nový AlphaFind CSV se vloží do:
#       data/raw/alphafind/inbox/
#
# - při prvním běhu se CSV automaticky přesune do nového:
#       data/raw/alphafind/projectN/
#
# - výstupy jdou do:
#       data/processed/projectN/
#
# - sdílené věci (NCBI taxonomy, taxon name cache, české názvy)
#   zůstávají mimo projectN a používají se pro všechny projekty.
#
# Běžný nový běh:
#       python tree.py
#
# Pokud je v inboxu právě jeden CSV soubor, vytvoří se automaticky
# další projectN.
#
# Opakovaný běh aktuálního projektu:
#       python tree.py --current
#
# Opakovaný běh konkrétního projektu:
#       python tree.py project2
#
# Nový projekt s hezkým názvem:
#       python tree.py --new --name "1FAS – fasciculin"


parser = argparse.ArgumentParser(
    description=(
        "Vytvoří taxonomické výstupy pro AlphaFind projekt. "
        "Technické projekty se jmenují project1, project2, ..."
    )
)

parser.add_argument(
    "project",
    nargs="?",
    help=(
        "Existující technické ID projektu, např. project2. "
        "Pokud není zadáno a v inboxu je CSV, vytvoří se nový projectN."
    )
)

parser.add_argument(
    "--new",
    action="store_true",
    help="Vynutí vytvoření nového projectN z CSV v inboxu."
)

parser.add_argument(
    "--current",
    action="store_true",
    help="Použije projekt uložený v data/current_project.txt."
)

parser.add_argument(
    "--name",
    help=(
        "Lidský název projektu ukládaný do project.json, "
        'např. "1FAS – fasciculin".'
    )
)

parser.add_argument(
    "--download",
    metavar="UNIPROT_ID",
    help=(
        "Stáhne výsledky přímo z AlphaFind API, uloží filtrované CSV "
        "do inboxu a vytvoří z něj nový projekt."
    )
)

parser.add_argument(
    "--index",
    default="chains",
    choices=["chains", "chains_90", "chains_80", "chains_70", "domains"],
    help="AlphaFind index pro --download (výchozí: chains)."
)

parser.add_argument(
    "--k",
    type=int,
    default=5000,
    help=(
        "Počet kNN kandidátů požadovaných z AlphaFind API "
        "(výchozí: 5000). Cut-off se aplikuje až na tyto kandidáty."
    )
)

parser.add_argument(
    "--tm-min",
    type=float,
    default=0.5,
    help="Minimální TM-score pro ponechání hitu (výchozí: 0.5)."
)

parser.add_argument(
    "--aligned-min",
    type=float,
    default=90.0,
    help=(
        "Minimální podíl aligned residues v procentech pro ponechání hitu "
        "(výchozí: 90)."
    )
)

parser.add_argument(
    "--no-hit-filter",
    action="store_true",
    help="Neaplikuje filtr TM-score / aligned residues na vstupní data."
)

parser.add_argument(
    "--csv-columns",
    choices=["auto", "correct", "swapped"],
    default="auto",
    help=(
        "Interpretace sloupců ručně exportovaného CSV. 'swapped' prohodí "
        "Aligned Residues a Sequence Identity před filtrováním a exportem. "
        "API data jsou vždy interpretována správně."
    )
)

parser.add_argument(
    "--poll-seconds",
    type=float,
    default=10.0,
    help="Interval kontroly výpočtu AlphaFind API (výchozí: 10 s)."
)

parser.add_argument(
    "--api-timeout-minutes",
    type=float,
    default=120.0,
    help="Nejdelší čekání na AlphaFind API (výchozí: 120 min)."
)

args = parser.parse_args()


# ------------------------------------------------------------
# ZÁKLADNÍ CESTY
# ------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"

ALPHAFIND_ROOT_DIR = RAW_DIR / "alphafind"
ALPHAFIND_INBOX_DIR = ALPHAFIND_ROOT_DIR / "inbox"

PROCESSED_ROOT_DIR = DATA_DIR / "processed"

# Stabilní index všech projektů pro viewer.
# Viewer díky němu nemusí procházet adresáře a hádat,
# které projectN složky existují.
PROJECT_INDEX_FILE = (
    PROCESSED_ROOT_DIR / "projects.json"
)

# Ukazatel na projekt, který byl vytvořen / použit naposledy.
CURRENT_PROJECT_FILE = DATA_DIR / "current_project.txt"


# ------------------------------------------------------------
# SDÍLENÉ CESTY NAPŘÍČ VŠEMI PROJEKTY
# ------------------------------------------------------------

# ČESKÉ NÁZVY TAXONŮ
TAXON_NAMES_CS_DIR = RAW_DIR / "taxon_names_cs"
TAXON_NAMES_CS_FILE = TAXON_NAMES_CS_DIR / "taxon_names_cs.tsv"

# CACHE NÁZVŮ TAXONŮ
# tree.py ji přímo nepotřebuje, ale cesta je zde záměrně popsána:
# cache patří mezi SDÍLENÁ data a nikdy se nepřesouvá do projectN.
TAXON_NAMES_CACHE_DIR = RAW_DIR / "taxon_names"
TAXON_NAMES_CACHE_FILE = TAXON_NAMES_CACHE_DIR / "taxon_names_cache.tsv"

# NCBI TAXONOMY
NCBI_TAXONOMY_DIR = RAW_DIR / "ncbi_taxonomy"

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDUMP_ARCHIVE = NCBI_TAXONOMY_DIR / "taxdump.tar.gz"

NODES_FILE = NCBI_TAXONOMY_DIR / "nodes.dmp"
NAMES_FILE = NCBI_TAXONOMY_DIR / "names.dmp"


# ------------------------------------------------------------
# ALPHAFIND API DOWNLOAD
# ------------------------------------------------------------

ALPHAFIND_API_URL = "https://alphafind.ics.muni.cz/api"
API_PAGE_SIZE = 100


def fraction_to_percent(value):
    """AlphaFind API ukládá aligned/identity jako 0–1; CSV používá %."""

    if value is None:
        return None

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None

    return numeric * 100.0 if abs(numeric) <= 1.0 else numeric


def api_result_to_row(result):
    """Převede jeden správně pojmenovaný API záznam na stabilní CSV řádek."""

    target_id = result.get("target_id") or result.get("id")

    # Chain výsledky mohou mít technický suffix _0; UniProt accession ne.
    protein_id = re.sub(r"_\d+$", "", str(target_id or ""))

    return {
        "Protein ID": protein_id,
        "Protein Name": result.get("protein_name"),
        "Gene Name": result.get("gene_name"),
        "Organism": result.get("organism"),
        "Tax ID": result.get("tax_id"),
        "TM-Score": result.get("tm_score_query"),
        "Target TM-Score": result.get("tm_score_target"),
        "RMSD": result.get("rmsd"),
        "Aligned Residues (%)": fraction_to_percent(
            result.get("aligned_residues")
        ),
        "Sequence Identity (%)": fraction_to_percent(
            result.get("sequential_identity")
        ),
        "kNN Score": result.get("score"),
        "Average pLDDT": result.get("avg_plddt")
    }


def fetch_result_page(session, query_id, index_name, page):
    """Načte jednu stránku; během výpočtu ji opakovaně kontroluje.

    Krátkodobé chyby serveru, spojení nebo odpovědi, která není validní JSON,
    nezastaví celou pipeline. Stránka se po pauze zkusí znovu až do vypršení
    celkového API timeoutu.
    """

    url = (
        f"{ALPHAFIND_API_URL}/search/{query_id}/"
        f"{index_name}/results"
    )

    deadline = time.monotonic() + args.api_timeout_minutes * 60.0

    while True:
        try:
            response = session.get(
                url,
                params={
                    "page": page,
                    "page_size": API_PAGE_SIZE,
                    "sort_by": "knn",
                    "sort_order": "desc"
                },
                timeout=120
            )

            response.raise_for_status()

            try:
                payload = response.json()
            except ValueError:
                preview = (response.text or "")[:200].replace("\n", " ")

                print(
                    f"  stránka {page}: server nevrátil platný JSON "
                    f"(HTTP {response.status_code}; "
                    f"odpověď: {preview!r}). "
                    f"Zkusím znovu za {args.poll_seconds:g} s."
                )

                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"AlphaFind stránka {page} opakovaně vracela "
                        "neplatnou odpověď."
                    )

                time.sleep(max(1.0, args.poll_seconds))
                continue

        except requests.RequestException as error:
            print(
                f"  stránka {page}: chyba spojení/API: {error}. "
                f"Zkusím znovu za {args.poll_seconds:g} s."
            )

            if time.monotonic() >= deadline:
                raise TimeoutError(
                    f"AlphaFind stránku {page} se nepodařilo stáhnout "
                    f"do {args.api_timeout_minutes:g} minut."
                ) from error

            time.sleep(max(1.0, args.poll_seconds))
            continue

        status = str(payload.get("status", "")).lower()
        results = payload.get("results") or []

        # U hotové stránky jsou všechny TM-score skutečně dopočítané.
        page_complete = (
            status == "completed"
            and all(
                result.get("tm_score_query") is not None
                and result.get("aligned_residues") is not None
                for result in results
            )
        )

        if page_complete:
            return payload

        if status in {"failed", "error"}:
            raise RuntimeError(
                f"AlphaFind výpočet selhal na stránce {page}: {payload}"
            )

        if time.monotonic() >= deadline:
            raise TimeoutError(
                f"AlphaFind výpočet nebyl hotový do "
                f"{args.api_timeout_minutes:g} minut (stránka {page})."
            )

        print(
            f"  stránka {page}: stav {status or 'čekám'}, "
            f"další kontrola za {args.poll_seconds:g} s"
        )
        time.sleep(max(1.0, args.poll_seconds))


def download_alphafind_to_inbox():
    """Spustí asynchronní hledání, stáhne všechny stránky a uloží CSV."""

    if not args.download:
        return None

    if args.project or args.current:
        parser.error("--download nelze kombinovat s projectN ani --current.")

    if args.k < 1:
        parser.error("--k musí být kladné celé číslo.")

    if args.tm_min < 0 or args.aligned_min < 0:
        parser.error("Cut-off hodnoty nesmějí být záporné.")

    ALPHAFIND_INBOX_DIR.mkdir(parents=True, exist_ok=True)

    existing = sorted(ALPHAFIND_INBOX_DIR.glob("*.csv"))
    if existing:
        raise RuntimeError(
            "Pro --download musí být AlphaFind inbox prázdný. "
            f"Teď obsahuje: {[path.name for path in existing]}"
        )

    query = args.download.strip()
    if not query:
        parser.error("--download vyžaduje neprázdné UniProt/PDB ID.")

    request_body = {
        "query": query,
        "index": [args.index],
        "options": {
            "k": args.k,
            "size": args.k
        }
    }

    print()
    print("Spouštím AlphaFind API hledání:")
    print(f"  query: {query}")
    print(f"  index: {args.index}")
    print(f"  kNN kandidátů: {args.k}")

    with requests.Session() as session:
        response = session.post(
            f"{ALPHAFIND_API_URL}/search",
            json=request_body,
            timeout=120
        )
        response.raise_for_status()
        search_payload = response.json()

        query_id = search_payload.get("id")
        if not query_id:
            raise RuntimeError(
                f"AlphaFind nevrátil query ID: {search_payload}"
            )

        print(f"  query_id: {query_id}")

        first_page = fetch_result_page(
            session,
            query_id,
            args.index,
            1
        )

        total_results = int(first_page.get("total_results") or 0)
        total_pages = int(
            first_page.get("total_pages")
            or math.ceil(total_results / API_PAGE_SIZE)
        )

        results = list(first_page.get("results") or [])

        for page in range(2, total_pages + 1):
            print(f"Stahuji AlphaFind stránku {page}/{total_pages}…")
            payload = fetch_result_page(
                session,
                query_id,
                args.index,
                page
            )
            results.extend(payload.get("results") or [])

    rows = [api_result_to_row(result) for result in results]
    downloaded_df = pd.DataFrame(rows)

    before_filter = len(downloaded_df)

    if not args.no_hit_filter:
        keep = (
            downloaded_df["TM-Score"].ge(args.tm_min)
            | downloaded_df["Aligned Residues (%)"].ge(args.aligned_min)
        )
        downloaded_df = downloaded_df.loc[keep].copy()

    safe_query = re.sub(r"[^A-Za-z0-9_.-]+", "_", query)
    filename = (
        f"alphafind_{safe_query}_{args.index}_k{args.k}_"
        f"tm{args.tm_min:g}_or_aligned{args.aligned_min:g}.csv"
    )
    output_path = ALPHAFIND_INBOX_DIR / filename

    downloaded_df.to_csv(
        output_path,
        sep=";",
        index=False
    )

    print()
    print(
        f"AlphaFind API: staženo {before_filter} kandidátů, "
        f"ponecháno {len(downloaded_df)} hitů."
    )
    print(
        f"Podmínka: TM-score >= {args.tm_min:g} NEBO "
        f"aligned residues >= {args.aligned_min:g} %."
    )
    print(f"CSV uloženo do: {output_path}")

    if total_results >= args.k:
        print()
        print(
            "UPOZORNĚNÍ: API vrátilo celý požadovaný počet k kandidátů. "
            "Mohou existovat další hity mimo kNN kandidátní množinu. "
            "Pro širší hledání spusť vyšší --k."
        )

    return {
        "query": query,
        "query_id": query_id,
        "index": args.index,
        "k": args.k,
        "downloaded_count": before_filter,
        "kept_count": len(downloaded_df),
        "tm_min": args.tm_min,
        "aligned_min": args.aligned_min,
        "column_mapping": "correct",
        "downloaded_at": datetime.now().isoformat()
    }


# ------------------------------------------------------------
# POMOCNÉ FUNKCE PRO PROJECT1, PROJECT2, ...
# ------------------------------------------------------------

PROJECT_RE = re.compile(r"^project(\d+)$", re.IGNORECASE)


def validate_project_id(project_id):
    """
    Povolené technické ID:
        project1
        project2
        project123

    Název složky je interní identifikátor a neměl by se ručně přejmenovávat.
    """

    if not PROJECT_RE.fullmatch(project_id):
        raise ValueError(
            f"Neplatné project ID: {project_id!r}. "
            "Použij tvar project1, project2, project3, ..."
        )

    return project_id.lower()


def project_number(project_id):
    match = PROJECT_RE.fullmatch(project_id)
    return int(match.group(1)) if match else None


def existing_project_ids():
    """
    Najde technické projectN složky jak v raw/alphafind,
    tak v processed/. Díky tomu nehrozí opětovné použití čísla,
    i kdyby jedna strana chyběla.
    """

    ids = set()

    for root in [ALPHAFIND_ROOT_DIR, PROCESSED_ROOT_DIR]:

        if not root.exists():
            continue

        for path in root.iterdir():

            if path.is_dir() and PROJECT_RE.fullmatch(path.name):
                ids.add(path.name.lower())

    return ids


def next_project_id():
    """
    Vrátí první další číselné ID:
        project1, project2, ...
    """

    used_numbers = [
        project_number(project_id)
        for project_id in existing_project_ids()
    ]

    used_numbers = [
        number
        for number in used_numbers
        if number is not None
    ]

    next_number = max(used_numbers, default=0) + 1

    return f"project{next_number}"


def read_current_project():
    """
    Přečte data/current_project.txt.
    """

    if not CURRENT_PROJECT_FILE.exists():
        return None

    value = CURRENT_PROJECT_FILE.read_text(
        encoding="utf-8"
    ).strip()

    if not value:
        return None

    return validate_project_id(value)


def write_current_project(project_id):
    """
    Zapamatuje technické ID projektu pro druhý průchod tree.py
    a pro další skripty v pipeline.
    """

    CURRENT_PROJECT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    CURRENT_PROJECT_FILE.write_text(
        project_id + "\n",
        encoding="utf-8"
    )


def inbox_csv_files():
    """
    Vrátí CSV soubory čekající v inboxu.
    """

    ALPHAFIND_INBOX_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    return sorted(
        ALPHAFIND_INBOX_DIR.glob("*.csv")
    )


def update_project_index(project_metadata):
    """
    Přidá nebo aktualizuje jeden projekt v:
        data/processed/projects.json

    Index je stabilní seznam, který používá viewer.

    Příklad:
    [
      {
        "project_id": "project1",
        "name": "Bet v 1",
        "input_csv": "...csv",
        "created_at": "...",
        "tree_json": "alphafind_tree.json"
      }
    ]
    """

    PROCESSED_ROOT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    projects = []

    if PROJECT_INDEX_FILE.exists():

        try:

            with open(
                PROJECT_INDEX_FILE,
                "r",
                encoding="utf-8"
            ) as f:

                loaded = json.load(f)

            if isinstance(loaded, list):
                projects = loaded

        except (
            json.JSONDecodeError,
            OSError
        ):

            print(
                "POZOR: projects.json se nepodařilo načíst. "
                "Vytvořím ho znovu."
            )

            projects = []


    # Do indexu ukládáme jen pole, která viewer potřebuje.
    entry = {
        "project_id":
            project_metadata.get(
                "project_id"
            ),

        "name":
            project_metadata.get(
                "name"
            )
            or project_metadata.get(
                "project_id"
            ),

        "input_csv":
            project_metadata.get(
                "input_csv"
            ),

        "created_at":
            project_metadata.get(
                "created_at"
            ),

        "tree_json":
            project_metadata.get(
                "tree_json",
                "alphafind_tree.json"
            )
    }


    # Starý záznam stejného project_id nahradíme.
    projects = [
        project
        for project in projects
        if project.get(
            "project_id"
        ) != entry[
            "project_id"
        ]
    ]

    projects.append(
        entry
    )


    def sort_key(project):

        project_id = str(
            project.get(
                "project_id",
                ""
            )
        )

        match = PROJECT_RE.fullmatch(
            project_id
        )

        if match:
            return int(
                match.group(1)
            )

        return 10**12


    projects.sort(
        key=sort_key
    )


    temp_file = (
        PROJECT_INDEX_FILE
        .with_suffix(".json.tmp")
    )


    with open(
        temp_file,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            projects,
            f,
            ensure_ascii=False,
            indent=2
        )


    temp_file.replace(
        PROJECT_INDEX_FILE
    )


def create_new_project(human_name=None):
    """
    Vytvoří další projectN.

    Očekává právě jeden CSV v:
        data/raw/alphafind/inbox/

    CSV se PŘESUNE do:
        data/raw/alphafind/projectN/

    Tím se inbox automaticky uvolní pro další AlphaFind export.
    """

    csv_files = inbox_csv_files()

    if len(csv_files) == 0:
        raise FileNotFoundError(
            "Chci vytvořit nový projekt, ale v inboxu není žádný CSV.\n"
            f"Vlož právě jeden AlphaFind CSV do:\n"
            f"  {ALPHAFIND_INBOX_DIR}"
        )

    if len(csv_files) > 1:
        raise RuntimeError(
            "V inboxu je více CSV souborů. "
            "Pro založení jednoho projektu tam nech právě jeden.\n"
            f"Inbox: {ALPHAFIND_INBOX_DIR}"
        )

    project_id = next_project_id()

    raw_project_dir = (
        ALPHAFIND_ROOT_DIR / project_id
    )

    processed_project_dir = (
        PROCESSED_ROOT_DIR / project_id
    )

    raw_project_dir.mkdir(
        parents=True,
        exist_ok=False
    )

    processed_project_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    source_csv = csv_files[0]

    destination_csv = (
        raw_project_dir / source_csv.name
    )

    shutil.move(
        str(source_csv),
        str(destination_csv)
    )

    project_metadata = {
        "project_id": project_id,

        # Toto je uživatelský název.
        # Můžeš ho kdykoliv změnit bez přejmenování složky.
        "name": (
            human_name.strip()
            if human_name and human_name.strip()
            else project_id
        ),

        "input_csv": source_csv.name,
        "created_at": datetime.now().isoformat(),

        # Relativní názvy hlavních výstupů.
        "tree_json": "alphafind_tree.json"
    }

    project_json_file = (
        processed_project_dir / "project.json"
    )

    with open(
        project_json_file,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            project_metadata,
            f,
            ensure_ascii=False,
            indent=2
        )

    update_project_index(
        project_metadata
    )

    write_current_project(project_id)

    print()
    print("Vytvořen nový projekt:")
    print(f"  ID:        {project_id}")
    print(f"  Název:     {project_metadata['name']}")
    print(f"  Raw CSV:   {destination_csv}")
    print(f"  Metadata:  {project_json_file}")

    return project_id


def load_or_create_project():
    """
    Rozhodne, který projekt použít.

    Priorita:
    1) explicitní projectN
    2) --current
    3) --new
    4) bez parametrů:
       - pokud je v inboxu jeden CSV -> nový projectN
       - pokud je inbox prázdný -> current_project.txt

    Tohle je záměrně vhodné i pro současný TreeRun.sh,
    který tree.py volá dvakrát:
       první průchod vytvoří projekt,
       druhý průchod už najde prázdný inbox a použije current project.
    """

    if args.project and args.current:
        parser.error(
            "Použij buď explicitní projectN, nebo --current, ne obojí."
        )

    if args.project and args.new:
        parser.error(
            "Použij buď explicitní projectN, nebo --new, ne obojí."
        )

    if args.current and args.new:
        parser.error(
            "Použij buď --current, nebo --new, ne obojí."
        )

    # 1) explicitní projectN
    if args.project:

        project_id = validate_project_id(
            args.project
        )

        write_current_project(project_id)

        return project_id

    # 2) aktuální projekt
    if args.current:

        project_id = read_current_project()

        if project_id is None:
            raise FileNotFoundError(
                "Soubor current_project.txt zatím neukazuje na žádný projekt."
            )

        return project_id

    # 3) explicitně nový projekt
    if args.new:

        return create_new_project(
            args.name
        )

    # 4) pohodlný automatický režim
    csv_files = inbox_csv_files()

    if len(csv_files) == 1:

        return create_new_project(
            args.name
        )

    if len(csv_files) > 1:

        raise RuntimeError(
            "V inboxu je více CSV souborů. "
            "Nech tam právě jeden."
        )

    # Inbox je prázdný -> znovu použij poslední projekt.
    project_id = read_current_project()

    if project_id is None:
        raise FileNotFoundError(
            "V inboxu není CSV a zároveň neexistuje current_project.txt.\n"
            f"Pro nový projekt vlož CSV do:\n"
            f"  {ALPHAFIND_INBOX_DIR}"
        )

    return project_id


# ------------------------------------------------------------
# URČIT AKTUÁLNÍ PROJEKT
# ------------------------------------------------------------

api_download_metadata = download_alphafind_to_inbox()

PROJECT_ID = load_or_create_project()

ALPHAFIND_DIR = (
    ALPHAFIND_ROOT_DIR / PROJECT_ID
)

PROCESSED_DIR = (
    PROCESSED_ROOT_DIR / PROJECT_ID
)

PROJECT_JSON_FILE = (
    PROCESSED_DIR / "project.json"
)

if not ALPHAFIND_DIR.exists():
    raise FileNotFoundError(
        f"Raw složka projektu neexistuje:\n"
        f"  {ALPHAFIND_DIR}"
    )

PROCESSED_DIR.mkdir(
    parents=True,
    exist_ok=True
)

TAXON_NAMES_CS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

TAXON_NAMES_CACHE_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ------------------------------------------------------------
# NAČÍST / DOPLNIT PROJECT.JSON
# ------------------------------------------------------------

if PROJECT_JSON_FILE.exists():

    with open(
        PROJECT_JSON_FILE,
        "r",
        encoding="utf-8"
    ) as f:

        project_metadata = json.load(f)

else:

    project_metadata = {
        "project_id": PROJECT_ID,
        "name": PROJECT_ID,
        "created_at": datetime.now().isoformat(),
        "tree_json": "alphafind_tree.json"
    }


# --name může změnit lidský název i při opakovaném běhu projektu.
if args.name and args.name.strip():

    project_metadata["name"] = (
        args.name.strip()
    )


project_metadata["project_id"] = PROJECT_ID
project_metadata.setdefault(
    "tree_json",
    "alphafind_tree.json"
)

if api_download_metadata is not None:
    project_metadata["alphafind_api"] = api_download_metadata

with open(
    PROJECT_JSON_FILE,
    "w",
    encoding="utf-8"
) as f:

    json.dump(
        project_metadata,
        f,
        ensure_ascii=False,
        indent=2
    )


update_project_index(
    project_metadata
)


write_current_project(PROJECT_ID)


print()
print("=" * 70)
print(f"PROJECT ID:   {PROJECT_ID}")
print(
    f"PROJECT NAME: "
    f"{project_metadata.get('name', PROJECT_ID)}"
)
print(f"AlphaFind:    {ALPHAFIND_DIR}")
print(f"Výstupy:      {PROCESSED_DIR}")
print(f"project.json: {PROJECT_JSON_FILE}")
print("=" * 70)


# ============================================================
# 2. INFORMACE O AKTUÁLNÍM OPEN TREE OF LIFE
# ============================================================

ABOUT_URL = "https://api.opentreeoflife.org/v3/tree_of_life/about"

response = requests.post(
    ABOUT_URL,
    json={"include_source_list": False},
    timeout=30,
)

response.raise_for_status()

info = response.json()

metadata = {
    "project_id": PROJECT_ID,
    "project_name": project_metadata.get("name", PROJECT_ID),
    "source": "Open Tree of Life",
    "downloaded_at": datetime.now().isoformat(),
    "synth_id": info["synth_id"],
    "taxonomy_version": info["taxonomy_version"],
    "date_created": info["date_created"],
    "root_ott_id": info["root"]["taxon"]["ott_id"],
    "root_name": info["root"]["taxon"]["name"],
    "num_tips": info["root"]["num_tips"],
}

METADATA_FILE = PROCESSED_DIR / "opentree_metadata.json"

with open(METADATA_FILE, "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print("OpenTree metadata:")
print(json.dumps(metadata, indent=2))

print()
print(f"Metadata saved to: {METADATA_FILE}")


# ============================================================
# 3. NAJÍT CSV EXPORT Z ALPHAFINDU
# ============================================================

csv_files = list(ALPHAFIND_DIR.glob("*.csv"))
#Ve složce ALPHAFIND_DIR najdi všechny soubory, které končí na .csv, a ulož je do seznamu csv_files.(* = libovolný název)
    #pomocí metody .glob (z objektu Path) hledá ve složce ALPHAFIND_DIR podle vzoru: *.csv, kde * znamená „cokoli před .csv“.
    #je v podstatě Pythonová obdoba bashového: ls *.csv akorát že python s tím pak může dál pracovat

if len(csv_files) == 0:
    raise FileNotFoundError(
        f"Ve složce datasetu {ALPHAFIND_DIR} není žádný CSV soubor."
    )

if len(csv_files) > 1:
    raise RuntimeError(
        f"Ve složce datasetu {ALPHAFIND_DIR} je více CSV souborů. "
        "Chcu jen jeden, kurňa, vymaž ty ostatní, je v tom pak bordel!."
    )

INPUT_FILE = csv_files[0]
# v Pythonu se položky v seznamu počítají od nuly! Vezme tedy první csv ze seznamu a tu uloží do INPUT_FILE

# Uložit název vstupního CSV také do projektových metadat.
project_metadata["input_csv"] = INPUT_FILE.name

with open(
    PROJECT_JSON_FILE,
    "w",
    encoding="utf-8"
) as f:

    json.dump(
        project_metadata,
        f,
        ensure_ascii=False,
        indent=2
    )


update_project_index(
    project_metadata
)


print()
print(f"Načítám AlphaFind CSV: {INPUT_FILE.name}")


# ============================================================
# 4. NAČÍST ALPHAFIND CSV DO PANDAS
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep=None,
    engine="python"
)
# pd je zkratka pro pandas; říká to "pandase, načti CSV soubor(již uložený pod proměnnou INPUT_FILE) a udělej z něj tabulku."
# sep=";" = sloupce jsou v souboru oddělené středníkem; když by byly jinak, tak napsat třeba čárku
# výsledek se uloží do df (běžná zkratka pro DataFrame); DataFrame je prostě pandasová tabulka.
print()
print("Prvních 5 řádků AlphaFind tabulky:")
print(df.head())
# () vypíše prázdný řádek
# napiš že tam bude "prvních 5 řádků..." a skutečně je zobraz (df.head ()); head() bez čísla defaultně ukazuje 5 řádků.
#   Můžu změnit třeba na (df.head(10))


# ============================================================
# 4B. NAJÍT SLOUPCE S ÚDAJI O STRUKTURNÍ PODOBNOSTI
# ============================================================

def normalize_column_name(column_name):
    """
    Převede různé zápisy názvu sloupce na podobný tvar.

    Např.:
    TM-score
    TM Score
    tm_score

    -> tm score
    """

    return " ".join(
        str(column_name)
        .strip()
        .lower()
        .replace("_", " ")
        .replace("-", " ")
        .split()
    )


def find_column(candidates, required=True):
    """
    Najde v AlphaFind CSV sloupec podle několika možných názvů.

    Pokud required=True a žádný kandidát neexistuje,
    skript skončí s chybou a vypíše skutečné názvy sloupců.
    """

    normalized_columns = {
        normalize_column_name(column): column
        for column in df.columns
    }

    for candidate in candidates:

        normalized_candidate = normalize_column_name(
            candidate
        )

        if normalized_candidate in normalized_columns:
            return normalized_columns[normalized_candidate]

    if required:
        raise KeyError(
            f"Nenašla jsem sloupec odpovídající: {candidates}\n"
            f"Sloupce v CSV jsou: {df.columns.tolist()}"
        )

    return None


# Základní strukturální údaje, které chceme určitě.
# V aktuálním AlphaFind exportu jsou sloupce:
#   TM-Score
#   Aligned Residues (%)

TM_SCORE_COL = find_column([
    "TM-Score",
    "TM-score",
    "TM score",
    "tm_score"
])


ALIGNED_PERCENT_COL = find_column([
    "Aligned Residues (%)",
    "Aligned percentage",
    "Aligned percent",
    "Aligned %",
    "aligned_percentage"
])


SEQUENCE_IDENTITY_COL = find_column([
    "Sequence Identity (%)",
    "Sequence identity",
    "Sequential Identity",
    "Sequential identity",
    "Sequence identity percent",
    "sequential_identity"
], required=False)


# Identifikace proteinu - kandidáti jsou schválně širší,
# protože různé AlphaFind exporty mohou používat jiný název.
PROTEIN_ID_COL = find_column([
    "Protein ID",
    "UniProt ID",
    "UniProt",
    "Object ID",
    "Entry",
    "Accession"
], required=False)


PROTEIN_NAME_COL = find_column([
    "Protein Name",
    "Protein name",
    "UniProt Description",
    "Description"
], required=False)


print()
print("Rozpoznané AlphaFind sloupce:")
print(f"TM-score: {TM_SCORE_COL}")
print(f"Aligned residues (%): {ALIGNED_PERCENT_COL}")
print(f"Sequence identity (%): {SEQUENCE_IDENTITY_COL}")
print(f"Protein ID: {PROTEIN_ID_COL}")
print(f"Protein name: {PROTEIN_NAME_COL}")


def clean_numeric_column(column):
    """
    Převede numerický sloupec z CSV na skutečná čísla.

    Umí odstranit znak % a převést desetinnou čárku na tečku.
    Neplatné hodnoty převede na NaN.
    """

    if column is None:
        return

    df[column] = pd.to_numeric(
        df[column]
        .astype(str)
        .str.replace("%", "", regex=False)
        .str.replace(",", ".", regex=False),
        errors="coerce"
    )


clean_numeric_column(TM_SCORE_COL)
clean_numeric_column(ALIGNED_PERCENT_COL)
clean_numeric_column(SEQUENCE_IDENTITY_COL)


# AlphaFind API má pole pojmenovaná správně. U ručního CSV exportu lze
# kvůli známé frontendové chybě explicitně zvolit --csv-columns swapped.
column_mapping = args.csv_columns

if project_metadata.get("alphafind_api"):
    column_mapping = "correct"
elif column_mapping == "auto":
    previous_mapping = (
        project_metadata
        .get("score_columns", {})
        .get("interpretation")
    )

    if previous_mapping in {"correct", "swapped"}:
        column_mapping = previous_mapping
        print(
            f"Používám dříve uloženou interpretaci CSV sloupců: "
            f"{column_mapping}."
        )
    else:
        column_mapping = "correct"
        print()
        print(
            "POZOR: --csv-columns auto u nového ručního CSV zachovává "
            "názvy sloupců. Pokud CSV pochází z chybné verze webu, "
            "použij při prvním běhu --csv-columns swapped."
        )

if column_mapping == "swapped":
    if SEQUENCE_IDENTITY_COL is None:
        raise KeyError(
            "Nelze použít --csv-columns swapped: CSV nemá sloupec "
            "Sequence Identity (%)."
        )

    aligned_values = df[ALIGNED_PERCENT_COL].copy()
    df[ALIGNED_PERCENT_COL] = df[SEQUENCE_IDENTITY_COL]
    df[SEQUENCE_IDENTITY_COL] = aligned_values

    print()
    print(
        "Pro filtrování a JSON export prohazuji Aligned Residues "
        "a Sequence Identity (--csv-columns swapped)."
    )


# Jednotné interní procentní měřítko. Některé exporty/API používají 0–1,
# jiné už 0–100. Hodnoty s absolutní velikostí <= 1 převedeme na procenta.
for percent_column in [ALIGNED_PERCENT_COL, SEQUENCE_IDENTITY_COL]:
    if percent_column is None:
        continue

    non_null = df[percent_column].dropna()
    if not non_null.empty and non_null.abs().max() <= 1.0:
        df[percent_column] = df[percent_column] * 100.0


# Ponecháváme hit, pokud splňuje ALESPOŇ jednu podmínku. To zachová i
# vícedoménový protein s nižším globálním TM-score, pokud se jedna doména
# zarovná přes vysoký podíl reziduí.
rows_before_filter = len(df)

if not args.no_hit_filter:
    keep_hit = (
        df[TM_SCORE_COL].ge(args.tm_min)
        | df[ALIGNED_PERCENT_COL].ge(args.aligned_min)
    )
    df = df.loc[keep_hit].copy()

    print()
    print(
        f"Filtr hitů: {rows_before_filter} -> {len(df)} řádků "
        f"(TM-score >= {args.tm_min:g} NEBO "
        f"aligned residues >= {args.aligned_min:g} %)."
    )

if df.empty:
    raise ValueError(
        "Po aplikaci filtru nezůstal žádný hit. "
        "Sniž --tm-min / --aligned-min nebo použij --no-hit-filter."
    )

project_metadata["score_columns"] = {
    "interpretation": column_mapping,
    "tm_min": None if args.no_hit_filter else args.tm_min,
    "aligned_min": None if args.no_hit_filter else args.aligned_min,
    "filter_operator": None if args.no_hit_filter else "or",
    "rows_before_filter": rows_before_filter,
    "rows_after_filter": len(df)
}

with open(PROJECT_JSON_FILE, "w", encoding="utf-8") as f:
    json.dump(project_metadata, f, ensure_ascii=False, indent=2)

update_project_index(project_metadata)


# ============================================================
# 5. VYTVOŘIT SOUHRN ORGANISMŮ + KVALITY HITŮ
# ============================================================

# Pro každý organismus chceme:
# - počet AlphaFind hitů
# - nejlepší / průměrný / medián TM-score
# - nejlepší / průměrný / medián podílu aligned residues (%)

aggregation = {

    "hit_count": (
        "Tax ID",
        "size"
    ),

    "tm_score_best": (
        TM_SCORE_COL,
        "max"
    ),

    "tm_score_mean": (
        TM_SCORE_COL,
        "mean"
    ),

    "tm_score_median": (
        TM_SCORE_COL,
        "median"
    ),

    "aligned_percent_best": (
        ALIGNED_PERCENT_COL,
        "max"
    ),

    "aligned_percent_mean": (
        ALIGNED_PERCENT_COL,
        "mean"
    ),

    "aligned_percent_median": (
        ALIGNED_PERCENT_COL,
        "median"
    )
}

if SEQUENCE_IDENTITY_COL is not None:
    aggregation.update({
        "sequence_identity_percent_best": (
            SEQUENCE_IDENTITY_COL,
            "max"
        ),
        "sequence_identity_percent_mean": (
            SEQUENCE_IDENTITY_COL,
            "mean"
        ),
        "sequence_identity_percent_median": (
            SEQUENCE_IDENTITY_COL,
            "median"
        )
    })


taxa = (
    df
    .groupby(
        ["Tax ID", "Organism"]
    )
    .agg(**aggregation)
    .reset_index()
    .sort_values(
        "hit_count",
        ascending=False
    )
)


# Zaokrouhlení, aby byl výstup čitelnější.

for column in [
    "tm_score_best",
    "tm_score_mean",
    "tm_score_median"
]:
    taxa[column] = taxa[column].round(3)


for column in [
    "aligned_percent_best",
    "aligned_percent_mean",
    "aligned_percent_median",
    "sequence_identity_percent_best",
    "sequence_identity_percent_mean",
    "sequence_identity_percent_median"
]:
    if column in taxa.columns:
        taxa[column] = taxa[column].round(1)


print()
print("Souhrn organismů a kvality hitů:")
print(taxa)

print()
print(
    f"Počet unikátních organismů: "
    f"{len(taxa)}"
)

# ============================================================
# 6. ULOŽIT ZPRACOVANOU TABULKU
# ============================================================

OUTPUT_FILE = PROCESSED_DIR / "alphafind_taxa.tsv"
# ulož do složky PROCESSED_DIR tady ten nový soubor nazvaný takhle a zároveň ho ulož jako proměnnou "OUTPUT_FILE"
#ještě je prázdný

taxa.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)
# vezmi pandas tabulku taxa a ulož ji do souboru OUTPUT_FILE.
# sloupce odděluj tabulátorem: \t znamená tabulátor (takže výsledný soubor je TSV)
#   Trochu matoucí je, že soubor ukládáme jako .tsv, ale používáme to_csv(). 
#   To je v pohodě — to_csv() umí ukládat tabulky s různými oddělovači.
#index=False znamená neukládej do souboru pandasový index řádků (0,1,2.. jsou interní indexy řádků co tam nechceme)

print()
print(f"Taxonomická tabulka uložena do: {OUTPUT_FILE}")
# vypíše cestu, kde byla tabulka uložena´

# ==========================================================
# ============================================================
# 7. STÁHNOUT A ROZBALIT NCBI TAXONOMII, POKUD JEŠTĚ NENÍ
# ============================================================

NCBI_TAXONOMY_DIR.mkdir(parents=True, exist_ok=True)

if not NODES_FILE.exists() or not NAMES_FILE.exists():

    print()
    print("NCBI taxonomy zatím není připravená.")

    # --------------------------------------------------------
    # stáhnout archiv, pokud ještě neexistuje
    # --------------------------------------------------------

    if not TAXDUMP_ARCHIVE.exists():

        print("Stahuji taxdump z NCBI...")

        response = requests.get(
            TAXDUMP_URL,
            stream=True,
            timeout=120
        )

        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))
        downloaded = 0

        with open(TAXDUMP_ARCHIVE, "wb") as f:

            for chunk in response.iter_content(chunk_size=1024 * 1024):

                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)

                    if total_size > 0:
                        print(
                            f"\rStaženo: {downloaded / 1024 / 1024:.1f} MB "
                            f"z {total_size / 1024 / 1024:.1f} MB",
                            end="",
                            flush=True
                        )

                    else:
                        print(
                            f"\rStaženo: {downloaded / 1024 / 1024:.1f} MB",
                            end="",
                            flush=True
                        )

        print()
        print("Taxdump stažen.")

    else:
        print("Archiv taxdump.tar.gz už existuje.")

    # --------------------------------------------------------
    # rozbalit archiv
    # --------------------------------------------------------

    print("Rozbaluji taxdump...")

    with tarfile.open(TAXDUMP_ARCHIVE, "r:gz") as tar:
        tar.extractall(NCBI_TAXONOMY_DIR)

    print("Taxdump rozbalen.")

else:

    print()
    print("NCBI taxonomy už je rozbalená, používám lokální kopii.")

    # ============================================================
# 8. NAČÍST NCBI TAXONOMII DO JEDNODUCHÝCH MAP
# ============================================================

parent = {}
rank = {}
name = {}
children = {}

# ------------------------------------------------------------
# nodes.dmp
# ------------------------------------------------------------
# Z každého řádku vezmeme:
# taxid, parent_taxid, rank
#
# a vytvoříme:
# parent[taxid] = parent_taxid
# rank[taxid] = rank
# children[parent_taxid] = seznam dětí
# ------------------------------------------------------------

with open(NODES_FILE, "r", encoding="utf-8") as f:

    for line in f:

        parts = line.split("|")

        taxid = int(parts[0].strip())
        parent_taxid = int(parts[1].strip())
        tax_rank = parts[2].strip()

        parent[taxid] = parent_taxid
        rank[taxid] = tax_rank

        if parent_taxid not in children:
            children[parent_taxid] = []

        children[parent_taxid].append(taxid)


# ------------------------------------------------------------
# names.dmp
# ------------------------------------------------------------
# Jeden taxon může mít více názvů:
# scientific name, synonym, common name...
#
# My chceme jen "scientific name".
# ------------------------------------------------------------

with open(NAMES_FILE, "r", encoding="utf-8") as f:

    for line in f:

        parts = line.split("|")

        taxid = int(parts[0].strip())
        tax_name = parts[1].strip()
        name_class = parts[3].strip()

        if name_class == "scientific name":
            name[taxid] = tax_name


print()
print(f"Načteno taxonů v parent mapě: {len(parent)}")
print(f"Načteno scientific names: {len(name)}")


# ============================================================
# 8B. NAČÍST ČESKÉ NÁZVY TAXONŮ
# ============================================================

# Soubor očekává dva sloupce:
#
# taxid    name_cs
# 3745     růžovité
# 3749     jabloň
# 3750     jabloň domácí
#
# Když soubor zatím neexistuje, skript ho vytvoří jen s hlavičkou
# a normálně pokračuje. To znamená, že češtinu můžeme doplňovat
# postupně, aniž bychom rozbili zbytek pipeline.

if not TAXON_NAMES_CS_FILE.exists():

    with open(
        TAXON_NAMES_CS_FILE,
        "w",
        encoding="utf-8"
    ) as f:

        f.write("taxid\tname_cs\n")

    print()
    print(
        "Soubor s českými názvy zatím neexistoval, "
        f"vytvořila jsem prázdnou šablonu: {TAXON_NAMES_CS_FILE}"
    )


cs_names_df = pd.read_csv(
    TAXON_NAMES_CS_FILE,
    sep="\t",
    dtype={
        "taxid": "Int64",
        "name_cs": "string"
    }
)


required_cs_columns = {
    "taxid",
    "name_cs"
}


if not required_cs_columns.issubset(
    cs_names_df.columns
):

    raise ValueError(
        f"Soubor {TAXON_NAMES_CS_FILE} musí obsahovat sloupce "
        "'taxid' a 'name_cs'."
    )


# Vyhodit prázdné řádky a vytvořit jednoduchou mapu:
#
# name_cs[3750] = "jabloň domácí"

cs_names_df = cs_names_df.dropna(
    subset=[
        "taxid",
        "name_cs"
    ]
)


cs_names_df["name_cs"] = (
    cs_names_df["name_cs"]
    .astype(str)
    .str.strip()
)


cs_names_df = cs_names_df[
    cs_names_df["name_cs"] != ""
]


name_cs = dict(
    zip(
        cs_names_df["taxid"].astype(int),
        cs_names_df["name_cs"]
    )
)


print()
print(
    f"Načteno českých názvů taxonů: "
    f"{len(name_cs)}"
)

# ============================================================
# ============================================================
# 9. FUNKCE PRO ZÍSKÁNÍ LINEAGE JEDNOHO TAXONU
# ============================================================

def get_lineage(taxid):

    lineage = []
    current = taxid

    while True:

        lineage.append(current)

        if parent[current] == current:
            break

        current = parent[current]

    return lineage


# ============================================================
# 10. SLOUČIT LINEAGE VŠECH EVIDENČNÍCH ORGANISMŮ
#
# AlphaFind hit-path zůstává oddělený od externích evidence vrstev.
# Taxonomická kostra se ale staví ze sjednocení:
#   AlphaFind + PANTHER + Pfam + InterPro family
# ============================================================

hit_path_taxa = set()

for _, row in taxa.iterrows():

    taxid = int(row["Tax ID"])

    if taxid not in parent:
        print(
            f"POZOR: Tax ID {taxid} "
            f"({row['Organism']}) není v NCBI taxonomy."
        )
        continue

    hit_path_taxa.update(get_lineage(taxid))


def load_external_taxa(tsv_path, source_label):
    """Načte přímé Tax ID z volitelné evidence vrstvy."""

    loaded = set()

    if not tsv_path.exists():
        print(f"{source_label}: evidence zatím není dostupná ({tsv_path.name}).")
        return loaded

    try:
        evidence_df = pd.read_csv(tsv_path, sep="\t")
    except Exception as error:
        print(f"POZOR: {source_label} evidence nelze načíst: {error}")
        return loaded

    if "taxid" not in evidence_df.columns:
        print(f"POZOR: {source_label} TSV nemá sloupec 'taxid': {tsv_path}")
        return loaded

    for value in evidence_df["taxid"].dropna():
        try:
            taxid = int(value)
        except (TypeError, ValueError):
            continue

        if taxid not in parent:
            print(
                f"POZOR: {source_label} Tax ID {taxid} "
                "není v lokální NCBI taxonomy."
            )
            continue

        loaded.add(taxid)

    print(f"{source_label}: načteno {len(loaded)} přímých taxonů.")
    return loaded


panther_direct_taxa = load_external_taxa(
    PROCESSED_DIR / "panther" / "panther_taxa.tsv",
    "PANTHER"
)

pfam_direct_taxa = load_external_taxa(
    PROCESSED_DIR / "pfam" / "pfam_taxa.tsv",
    "Pfam"
)

family_direct_taxa = load_external_taxa(
    PROCESSED_DIR / "family" / "family_taxa.tsv",
    "InterPro family"
)

# Kostra stromu = AlphaFind + externí evidence + všichni jejich předci.
skeleton_path_taxa = set(hit_path_taxa)

for taxid in sorted(
    panther_direct_taxa
    | pfam_direct_taxa
    | family_direct_taxa
):
    skeleton_path_taxa.update(
        get_lineage(taxid)
    )

print()
print(f"AlphaFind hit-path taxonů: {len(hit_path_taxa)}")
print(f"PANTHER přímých taxonů: {len(panther_direct_taxa)}")
print(f"Pfam přímých taxonů: {len(pfam_direct_taxa)}")
print(f"InterPro family přímých taxonů: {len(family_direct_taxa)}")
print(f"Taxonů v evidenční kostře včetně předků: {len(skeleton_path_taxa)}")

# ============================================================
# ============================================================
# 11. PŘIDAT SESTERSKÉ TAXONY JAKO KONTEXT
#     vytvoříme FULL a COMPACT verzi stromu
# ============================================================

from collections import Counter


# ------------------------------------------------------------
# 11A. FULL CONTEXT
# ------------------------------------------------------------

full_context_taxa = set()

for taxid in skeleton_path_taxa:

    for child_taxid in children.get(taxid, []):

        if child_taxid not in skeleton_path_taxa:
            full_context_taxa.add(child_taxid)


full_tree_taxa = skeleton_path_taxa | full_context_taxa


print()
print("FULL TREE:")
print(f"Evidenční kostra taxonů: {len(skeleton_path_taxa)}")
print(f"Kontextových sesterských taxonů: {len(full_context_taxa)}")
print(f"Celkem taxonů ve stromu: {len(full_tree_taxa)}")


# ------------------------------------------------------------
# DIAGNOSTIKA FULL CONTEXT PODLE RANKU
# ------------------------------------------------------------

full_context_ranks = Counter(
    rank.get(taxid, "UNKNOWN")
    for taxid in full_context_taxa
)

print()
print("FULL CONTEXT - taxony podle ranku:")

for tax_rank, count in full_context_ranks.most_common():
    print(f"{tax_rank}: {count}")


# ------------------------------------------------------------
# 11B. COMPACT CONTEXT
# ------------------------------------------------------------

LOW_LEVEL_RANKS = {
    "species",
    "subspecies",
    "varietas",
    "forma",
    "strain"
}

compact_context_taxa = set()

for taxid in skeleton_path_taxa:

    for child_taxid in children.get(taxid, []):

        if child_taxid in skeleton_path_taxa:
            continue

        child_rank = rank.get(child_taxid, "UNKNOWN")

        if child_rank in LOW_LEVEL_RANKS:
            continue

        compact_context_taxa.add(child_taxid)


compact_tree_taxa = skeleton_path_taxa | compact_context_taxa


print()
print("COMPACT TREE:")
print(f"Evidenční kostra taxonů: {len(skeleton_path_taxa)}")
print(f"Kontextových sesterských taxonů: {len(compact_context_taxa)}")
print(f"Celkem taxonů ve stromu: {len(compact_tree_taxa)}")


# ------------------------------------------------------------
# DIAGNOSTIKA COMPACT CONTEXT PODLE RANKU
# ------------------------------------------------------------

compact_context_ranks = Counter(
    rank.get(taxid, "UNKNOWN")
    for taxid in compact_context_taxa
)

print()
print("COMPACT CONTEXT - taxony podle ranku:")

for tax_rank, count in compact_context_ranks.most_common():
    print(f"{tax_rank}: {count}")


# ------------------------------------------------------------
# DIAGNOSTIKA:
# KTEŘÍ RODIČE PŘIDÁVAJÍ NEJVÍC KONTEXTOVÝCH DĚTÍ
# ------------------------------------------------------------

context_children_per_parent = []

for taxid in skeleton_path_taxa:

    full_context_children = [
        child_taxid
        for child_taxid in children.get(taxid, [])
        if child_taxid not in skeleton_path_taxa
    ]

    if full_context_children:

        context_children_per_parent.append(
            (
                len(full_context_children),
                taxid,
                name.get(taxid, "UNKNOWN"),
                rank.get(taxid, "UNKNOWN")
            )
        )


context_children_per_parent.sort(reverse=True)


print()
print("TOP 20 uzlů s největším počtem kontextových dětí:")

for count, taxid, tax_name, tax_rank in context_children_per_parent[:20]:

    print(
        count,
        taxid,
        tax_name,
        tax_rank
    )

# ============================================================
# ============================================================
# 12. VYTVOŘIT TABULKY UZLŮ PRO FULL A COMPACT STROM
# ============================================================

hit_count_by_taxid = dict(
    zip(
        taxa["Tax ID"].astype(int),
        taxa["hit_count"].astype(int)
    )
)


def make_tree_nodes_table(tree_taxa, context_taxa):

    rows = []

    for taxid in tree_taxa:

        if taxid in hit_count_by_taxid:
            node_type = "hit"

        elif taxid in context_taxa:
            node_type = "context"

        else:
            node_type = "path"

        rows.append(
            {
                "taxid": taxid,
                "parent_taxid": parent.get(taxid),
                "name": name.get(taxid, "UNKNOWN"),
                "name_cs": name_cs.get(taxid),
                "rank": rank.get(taxid, "UNKNOWN"),
                "node_type": node_type,
                "hit_count": hit_count_by_taxid.get(taxid, 0)
            }
        )

    return pd.DataFrame(rows)


full_tree_nodes = make_tree_nodes_table(
    full_tree_taxa,
    full_context_taxa
)

FULL_NODES_FILE = PROCESSED_DIR / "alphafind_tree_nodes_full.tsv"

full_tree_nodes.to_csv(
    FULL_NODES_FILE,
    sep="\t",
    index=False
)

print()
print(f"FULL tree nodes uloženy do: {FULL_NODES_FILE}")


compact_tree_nodes = make_tree_nodes_table(
    compact_tree_taxa,
    compact_context_taxa
)

COMPACT_NODES_FILE = PROCESSED_DIR / "alphafind_tree_nodes_compact.tsv"

compact_tree_nodes.to_csv(
    COMPACT_NODES_FILE,
    sep="\t",
    index=False
)

print(f"COMPACT tree nodes uloženy do: {COMPACT_NODES_FILE}")

print()
print("COMPACT TREE - počty podle node_type:")
print(compact_tree_nodes["node_type"].value_counts())

# ============================================================
# ============================================================
# 13. VYTVOŘIT NEWICK + LABELS + COLLAPSE + HIT COLORS
# ============================================================


# ------------------------------------------------------------
# 13A. BEZPEČNÉ ID UZLU
# ------------------------------------------------------------

def make_node_id(taxid):
    """
    Vytvoří jednoznačné ID uzlu vhodné pro Newick a iTOL.

    Příklad:
    NCBI_3750
    """

    return f"NCBI_{taxid}"


# ------------------------------------------------------------
# 13B. MAPA PARENT -> CHILDREN
# ------------------------------------------------------------

def make_tree_children(tree_taxa):
    """
    Z množiny taxonů vytvoří mapu:

    parent_taxid -> [child_taxid, child_taxid, ...]

    Použijeme jen vztahy uvnitř našeho vybraného stromu.
    """

    tree_children = {}

    for taxid in tree_taxa:

        parent_taxid = parent[taxid]

        # NCBI root má vztah:
        # 1 -> 1
        if taxid == parent_taxid:
            continue

        # Rodič musí být součástí našeho stromu.
        if parent_taxid not in tree_taxa:
            continue

        if parent_taxid not in tree_children:
            tree_children[parent_taxid] = []

        tree_children[parent_taxid].append(taxid)

    return tree_children


# ------------------------------------------------------------
# 13C. REKURZIVNÍ PŘEVOD DO NEWICKU
# ------------------------------------------------------------

def taxon_to_newick(taxid, tree_children):
    """
    Rekurzivně převede taxon a jeho potomky
    do Newick zápisu.

    V Newicku používáme pouze bezpečná ID:
    NCBI_...
    """

    node_id = make_node_id(taxid)

    child_taxa = tree_children.get(taxid, [])

    # Pokud uzel nemá děti, je to leaf.
    if len(child_taxa) == 0:
        return node_id

    # Děti seřadíme podle scientific name,
    # aby byl výstup při každém běhu stejný.
    child_taxa = sorted(
        child_taxa,
        key=lambda child: name.get(child, "")
    )

    child_newicks = []

    for child_taxid in child_taxa:

        child_newick = taxon_to_newick(
            child_taxid,
            tree_children
        )

        child_newicks.append(child_newick)

    children_text = ",".join(child_newicks)

    return f"({children_text}){node_id}"


def make_newick_tree(tree_taxa, root_taxid=1):
    """
    Vytvoří celý Newick strom
    a přidá závěrečný středník.
    """

    tree_children = make_tree_children(tree_taxa)

    if root_taxid not in tree_taxa:
        raise ValueError(
            f"Root Tax ID {root_taxid} není součástí stromu."
        )

    newick = taxon_to_newick(
        root_taxid,
        tree_children
    )

    return newick + ";"


# ------------------------------------------------------------
# 13D. FULL NEWICK
# ------------------------------------------------------------

full_newick = make_newick_tree(
    full_tree_taxa,
    root_taxid=1
)

FULL_NEWICK_FILE = (
    PROCESSED_DIR / "alphafind_tree_full.nwk"
)

with open(
    FULL_NEWICK_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(full_newick)

print()
print(
    f"FULL Newick uložen do: {FULL_NEWICK_FILE}"
)


# ------------------------------------------------------------
# 13E. COMPACT NEWICK
# ------------------------------------------------------------

compact_newick = make_newick_tree(
    compact_tree_taxa,
    root_taxid=1
)

COMPACT_NEWICK_FILE = (
    PROCESSED_DIR / "alphafind_tree_compact.nwk"
)

with open(
    COMPACT_NEWICK_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(compact_newick)

print(
    f"COMPACT Newick uložen do: {COMPACT_NEWICK_FILE}"
)


# ------------------------------------------------------------
# 13F. iTOL LABELS
# ------------------------------------------------------------

LABELS_FILE = (
    PROCESSED_DIR / "alphafind_tree_labels.txt"
)

with open(
    LABELS_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write("LABELS\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATA\n")

    for taxid in sorted(
        compact_tree_taxa,
        key=lambda x: name.get(x, "")
    ):

        node_id = make_node_id(taxid)
        tax_name = name.get(taxid, "UNKNOWN")

        f.write(
            f"{node_id}\t{tax_name}\n"
        )

print(
    f"iTOL LABELS uloženy do: {LABELS_FILE}"
)


# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# 13G. iTOL COLLAPSE - DEFAULT VIEW
# ------------------------------------------------------------

COLLAPSE_FILE = (
    PROCESSED_DIR / "alphafind_tree_collapse.txt"
)

VIRIDIPLANTAE_TAXID = 33090

collapse_taxa = [
    VIRIDIPLANTAE_TAXID
]

with open(
    COLLAPSE_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write("COLLAPSE\n")
    f.write("DATA\n")

    for taxid in collapse_taxa:
        f.write(
            f"{make_node_id(taxid)}\n"
        )

print(
    f"iTOL default COLLAPSE uložen do: "
    f"{COLLAPSE_FILE}"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
# 13H. iTOL TREE COLORS
#      - přímé hity červeně a silněji
#      - všechny uzly vedoucí k hitům také červeně
# ------------------------------------------------------------

HIT_COLORS_FILE = (
    PROCESSED_DIR / "alphafind_tree_hit_colors.txt"
)

with open(
    HIT_COLORS_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write("TREE_COLORS\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATA\n")

    for taxid in sorted(hit_path_taxa):

        node_id = make_node_id(taxid)

        # ----------------------------------------------------
        # Přímý AlphaFind hit
        # ----------------------------------------------------

        if taxid in hit_count_by_taxid:

            # silnější červená větev
            f.write(
                f"{node_id}\tbranch\t#FF0000\tnormal\t3\n"
            )

            # červený a tučný label
            f.write(
                f"{node_id}\tlabel\t#FF0000\tbold\n"
            )

        # ----------------------------------------------------
        # Taxon na cestě k hitu
        # ----------------------------------------------------

        else:

            # tenčí červená větev
            f.write(
                f"{node_id}\tbranch\t#FF0000\tnormal\t2\n"
            )

            # i název nadřazené skupiny bude červený
            f.write(
                f"{node_id}\tlabel\t#FF0000\tnormal\n"
            )


print(
    f"iTOL HIT COLORS uloženy do: "
    f"{HIT_COLORS_FILE}"
)
# ------------------------------------------------------------
# 13I. DIAGNOSTIKA
# ------------------------------------------------------------

print()
print("Velikost Newick stromů:")

print(
    f"FULL: {len(full_newick):,} znaků"
)

print(
    f"COMPACT: {len(compact_newick):,} znaků"
)

print()
print(
    f"Počet uzlů v COLLAPSE souboru: "
    f"{len(collapse_taxa)}"
)

print(
    f"Počet hitových taxonů obarvených červeně: "
    f"{len(hit_count_by_taxid)}"
)

print()
print("Začátek COMPACT Newicku:")

print(
    compact_newick[:500]
)

# ============================================================
# 14. EXPORT COMPACT STROMU DO JSONU PRO INTERAKTIVNÍ VIEWER
#     + TM-SCORE + ALIGNED RESIDUES (%) + DETAILY HITŮ
# ============================================================

TREE_JSON_FILE = (
    PROCESSED_DIR / "alphafind_tree.json"
)

TREE_ROOT_TAXID = 1


# ------------------------------------------------------------
# 14A. CHILDREN MAPA PRO COMPACT STROM
# ------------------------------------------------------------

json_tree_children = make_tree_children(
    compact_tree_taxa
)


# ------------------------------------------------------------
# 14B. SPOČÍTAT POČET HITŮ V CELÉ PODVĚTVI
# ------------------------------------------------------------

def get_subtree_hit_count(taxid):
    """
    Spočítá všechny AlphaFind hity v daném taxonu
    a ve všech jeho potomcích v compact stromu.
    """

    total = hit_count_by_taxid.get(
        taxid,
        0
    )

    for child_taxid in json_tree_children.get(
        taxid,
        []
    ):
        total += get_subtree_hit_count(
            child_taxid
        )

    return total


# ------------------------------------------------------------
# 14C. POMOCNÉ FUNKCE PRO BEZPEČNÝ EXPORT DO JSON
# ------------------------------------------------------------

def safe_float(value):
    """
    NaN -> None
    číslo -> float
    """

    if pd.isna(value):
        return None

    return float(value)


def safe_text(value):
    """
    NaN -> None
    text -> str
    """

    if pd.isna(value):
        return None

    return str(value)


# ------------------------------------------------------------
# 14D. DETAILY JEDNOTLIVÝCH ALPHAFIND HITŮ
# ------------------------------------------------------------

# Výsledkem bude např.:
#
# hits_by_taxid[3750] = [
#     {
#         "protein_id": "...",
#         "protein_name": "...",
#         "tm_score": 0.87,
#         "aligned_percent": 94.2
#     },
#     ...
# ]

hits_by_taxid = {}


for taxid, group in df.groupby("Tax ID"):

    taxid = int(taxid)

    hits = []


    for _, row in group.iterrows():

        hit = {

            "tm_score": safe_float(
                row[TM_SCORE_COL]
            ),

            "aligned_percent": safe_float(
                row[ALIGNED_PERCENT_COL]
            ),

            "sequence_identity_percent": (
                safe_float(row[SEQUENCE_IDENTITY_COL])
                if SEQUENCE_IDENTITY_COL is not None
                else None
            )
        }


        if PROTEIN_ID_COL is not None:

            hit["protein_id"] = safe_text(
                row[PROTEIN_ID_COL]
            )


        if PROTEIN_NAME_COL is not None:

            hit["protein_name"] = safe_text(
                row[PROTEIN_NAME_COL]
            )


        hits.append(hit)


    # Nejlepší TM-score bude ve vieweru nahoře.
    hits.sort(
        key=lambda hit: (
            hit["tm_score"]
            if hit["tm_score"] is not None
            else -1
        ),
        reverse=True
    )


    hits_by_taxid[taxid] = hits


# ------------------------------------------------------------
# 14E. SOUHRNNÉ SCORE PODLE TAX ID
# ------------------------------------------------------------

taxa_scores_by_taxid = {}


for _, row in taxa.iterrows():

    taxid = int(
        row["Tax ID"]
    )

    taxa_scores_by_taxid[taxid] = {

        "tm_score_best": safe_float(
            row["tm_score_best"]
        ),

        "tm_score_mean": safe_float(
            row["tm_score_mean"]
        ),

        "tm_score_median": safe_float(
            row["tm_score_median"]
        ),

        "aligned_percent_best": safe_float(
            row["aligned_percent_best"]
        ),

        "aligned_percent_mean": safe_float(
            row["aligned_percent_mean"]
        ),

        "aligned_percent_median": safe_float(
            row["aligned_percent_median"]
        ),

        "sequence_identity_percent_best": (
            safe_float(row["sequence_identity_percent_best"])
            if "sequence_identity_percent_best" in row.index
            else None
        ),

        "sequence_identity_percent_mean": (
            safe_float(row["sequence_identity_percent_mean"])
            if "sequence_identity_percent_mean" in row.index
            else None
        ),

        "sequence_identity_percent_median": (
            safe_float(row["sequence_identity_percent_median"])
            if "sequence_identity_percent_median" in row.index
            else None
        )
    }


# ------------------------------------------------------------
# 14F. REKURZIVNĚ POSTAVIT JSON UZEL
# ------------------------------------------------------------

def build_json_node(taxid):
    """
    Vytvoří jeden uzel JSON stromu včetně potomků.

    U přímého AlphaFind hitu přidá také:
    - souhrn TM-score
    - souhrn Aligned Residues (%)
    - seznam jednotlivých proteinových hitů
    """

    child_taxa = json_tree_children.get(
        taxid,
        []
    )


    child_taxa = sorted(
        child_taxa,
        key=lambda child: name.get(
            child,
            ""
        )
    )


    direct_hit_count = hit_count_by_taxid.get(
        taxid,
        0
    )


    contains_hit = (
        taxid in hit_path_taxa
    )


    if taxid in hit_count_by_taxid:

        node_type = "hit"

    elif taxid in hit_path_taxa:

        node_type = "path"

    else:

        node_type = "context"


    score_summary = taxa_scores_by_taxid.get(
        taxid,
        {}
    )


    node = {

        "taxid": taxid,

        "name": name.get(
            taxid,
            "UNKNOWN"
        ),

        # český název, pokud ho máme v taxon_names_cs.tsv;
        # jinak bude v JSONu null
        "name_cs": name_cs.get(
            taxid
        ),

        "rank": rank.get(
            taxid,
            "UNKNOWN"
        ),

        "node_type": node_type,


        # počet AlphaFind hitů přímo v tomto taxonu
        "direct_hit_count": direct_hit_count,


        # zda taxon leží na cestě k alespoň jednomu hitu
        "contains_hit": contains_hit,


        # počet hitů v celé podvětvi
        "subtree_hit_count": get_subtree_hit_count(
            taxid
        ),


        # ----------------------------------------------------
        # SOUHRNNÉ ÚDAJE O STRUKTURNÍ PODOBNOSTI
        # ----------------------------------------------------

        "tm_score_best": score_summary.get(
            "tm_score_best"
        ),

        "tm_score_mean": score_summary.get(
            "tm_score_mean"
        ),

        "tm_score_median": score_summary.get(
            "tm_score_median"
        ),

        "aligned_percent_best": score_summary.get(
            "aligned_percent_best"
        ),

        "aligned_percent_mean": score_summary.get(
            "aligned_percent_mean"
        ),

        "aligned_percent_median": score_summary.get(
            "aligned_percent_median"
        ),

        "sequence_identity_percent_best": score_summary.get(
            "sequence_identity_percent_best"
        ),

        "sequence_identity_percent_mean": score_summary.get(
            "sequence_identity_percent_mean"
        ),

        "sequence_identity_percent_median": score_summary.get(
            "sequence_identity_percent_median"
        ),


        # ----------------------------------------------------
        # JEDNOTLIVÉ ALPHAFIND HITY TOHOTO TAXONU
        # ----------------------------------------------------

        "hits": hits_by_taxid.get(
            taxid,
            []
        ),


        "has_children": len(
            child_taxa
        ) > 0,


        "children": [
            build_json_node(
                child_taxid
            )
            for child_taxid in child_taxa
        ]
    }


    return node


# ------------------------------------------------------------
# 14G. VYTVOŘIT CELÝ STROM
# ------------------------------------------------------------

if TREE_ROOT_TAXID not in compact_tree_taxa:

    raise ValueError(
        f"Root Tax ID {TREE_ROOT_TAXID} "
        f"není součástí compact stromu."
    )


tree_json = build_json_node(
    TREE_ROOT_TAXID
)


# ------------------------------------------------------------
# 14H. ULOŽIT JSON
# ------------------------------------------------------------

with open(
    TREE_JSON_FILE,
    "w",
    encoding="utf-8"
) as f:

    json.dump(
        tree_json,
        f,
        ensure_ascii=False,
        indent=2
    )


print()
print(
    f"JSON strom uložen do: "
    f"{TREE_JSON_FILE}"
)


# ------------------------------------------------------------
# 14I. DIAGNOSTIKA
# ------------------------------------------------------------

print(
    f"Root: "
    f"{tree_json['name']} "
    f"(Tax ID {tree_json['taxid']})"
)


print(
    f"Hitů v celém stromu: "
    f"{tree_json['subtree_hit_count']}"
)


print(
    f"Přímých dětí rootu: "
    f"{len(tree_json['children'])}"
)


print()
print("Kontrola detailů hitů:")


hit_taxa_with_details = [
    taxid
    for taxid, hits in hits_by_taxid.items()
    if len(hits) > 0
]


print(
    f"Taxonů se seznamem jednotlivých hitů: "
    f"{len(hit_taxa_with_details)}"
)


if hit_taxa_with_details:

    example_taxid = hit_taxa_with_details[0]

    print(
        f"Příklad: Tax ID {example_taxid}, "
        f"{len(hits_by_taxid[example_taxid])} hitů"
    )

    print("První hit:")

    print(
        hits_by_taxid[example_taxid][0]
    )



print()
print("Kontrola českých názvů:")
print(
    f"Českých názvů načteno celkem: "
    f"{len(name_cs)}"
)

example_cs_taxa = [
    taxid
    for taxid in sorted(compact_tree_taxa)
    if taxid in name_cs
][:10]

if example_cs_taxa:

    print("Příklady českých názvů ve stromu:")

    for taxid in example_cs_taxa:

        print(
            f"  {taxid}: "
            f"{name.get(taxid, 'UNKNOWN')} "
            f"-> {name_cs[taxid]}"
        )

else:

    print(
        "V aktuálním stromu zatím není žádný český název. "
        "To je v pořádku, pokud je taxon_names_cs.tsv zatím prázdný."
    )


# ============================================================
# 15. ZÁVĚREČNÝ PŘEHLED PROJEKTU
# ============================================================

print()
print("=" * 70)
print("HOTOVO")
print(f"Project ID:    {PROJECT_ID}")
print(
    f"Název:         "
    f"{project_metadata.get('name', PROJECT_ID)}"
)
print(f"Vstupní CSV:   {INPUT_FILE}")
print(f"Výstupní složka: {PROCESSED_DIR}")
print(f"Metadata:      {PROJECT_JSON_FILE}")
print("=" * 70)
