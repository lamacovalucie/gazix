#!/usr/bin/env python3

# ============================================================
# prepare_taxon_names_cs_cached.py
#
# Účel:
# - zjistí aktuální projectN z data/current_project.txt
# - vezme taxony z:
#       data/processed/projectN/alphafind_tree_nodes_compact.tsv
# - zkontroluje sdílenou TSV cache
# - na Wikidata se ptá jen na taxony, které v cache ještě nebyly zkontrolovány
# - do SDÍLENÉ cache ukládá jak nalezené názvy, tak i informaci "not_found"
# - z cache vytvoří SDÍLENÝ český soubor taxon_names_cs.tsv
#
# Důležité:
# - cache je perzistentní mezi všemi projectN
# - český soubor je také společný pro všechny projectN
# - pouze COMPACT_NODES_FILE se mění podle aktuálního projektu
# - při HTTP 429 používáme retry + delší čekání
# ============================================================

from pathlib import Path
import time

import pandas as pd
import requests


# ============================================================
# 1. CESTY
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_ROOT_DIR = DATA_DIR / "processed"

# ------------------------------------------------------------
# AKTUÁLNÍ PROJEKT
# ------------------------------------------------------------

CURRENT_PROJECT_FILE = (
    DATA_DIR / "current_project.txt"
)

if not CURRENT_PROJECT_FILE.exists():
    raise FileNotFoundError(
        "\nChybí soubor:\n"
        f"{CURRENT_PROJECT_FILE}\n\n"
        "Ten vytváří tree.py. "
        "Nejdřív spusť tree.py, aby bylo jasné, "
        "se kterým projectN se právě pracuje."
    )

PROJECT_ID = (
    CURRENT_PROJECT_FILE
    .read_text(encoding="utf-8")
    .strip()
)

if not PROJECT_ID:
    raise ValueError(
        f"{CURRENT_PROJECT_FILE} je prázdný."
    )

# Dataset-specific cesta:
PROCESSED_DIR = (
    PROCESSED_ROOT_DIR / PROJECT_ID
)

COMPACT_NODES_FILE = (
    PROCESSED_DIR
    / "alphafind_tree_nodes_compact.tsv"
)


# ------------------------------------------------------------
# SDÍLENÉ SOUBORY
# ------------------------------------------------------------

TAXON_NAMES_DIR = (
    RAW_DIR / "taxon_names"
)

TAXON_NAMES_CACHE_FILE = (
    TAXON_NAMES_DIR
    / "taxon_names_cache.tsv"
)

TAXON_NAMES_CS_DIR = (
    RAW_DIR / "taxon_names_cs"
)

TAXON_NAMES_CS_FILE = (
    TAXON_NAMES_CS_DIR
    / "taxon_names_cs.tsv"
)


print()
print("=" * 60)
print(f"Aktuální projekt: {PROJECT_ID}")
print(f"COMPACT strom:    {COMPACT_NODES_FILE}")
print(f"Cache:            {TAXON_NAMES_CACHE_FILE}")
print("=" * 60)


# ============================================================
# 2. NASTAVENÍ WIKIDAT
# ============================================================

WIKIDATA_SPARQL_URL = (
    "https://query.wikidata.org/sparql"
)

LANGUAGE = "cs"

# Menší dávky jsou šetrnější k veřejnému endpointu.
BATCH_SIZE = 50

# Pauza mezi úspěšnými dávkami.
PAUSE_SECONDS = 1.5

# Retry při 429 / dočasné síťové chybě.
MAX_RETRIES = 4

# Čekání se bude násobit 2**attempt:
# 5 s, 10 s, 20 s, 40 s
RETRY_BASE_SECONDS = 5

HEADERS = {
    "User-Agent": (
        "Lucka-Structures taxonomy viewer "
        "(research/educational use)"
    ),
    "Accept": "application/sparql-results+json"
}


# ============================================================
# 3. PŘIPRAVIT SLOŽKY A SOUBORY
# ============================================================

TAXON_NAMES_DIR.mkdir(
    parents=True,
    exist_ok=True
)

TAXON_NAMES_CS_DIR.mkdir(
    parents=True,
    exist_ok=True
)


if not TAXON_NAMES_CACHE_FILE.exists():

    pd.DataFrame(
        columns=[
            "taxid",
            "lang",
            "label",
            "source",
            "status"
        ]
    ).to_csv(
        TAXON_NAMES_CACHE_FILE,
        sep="\t",
        index=False
    )

    print(
        "Vytvořena prázdná cache: "
        f"{TAXON_NAMES_CACHE_FILE}"
    )


if not TAXON_NAMES_CS_FILE.exists():

    pd.DataFrame(
        columns=[
            "taxid",
            "name_cs"
        ]
    ).to_csv(
        TAXON_NAMES_CS_FILE,
        sep="\t",
        index=False
    )

    print(
        "Vytvořen český names soubor: "
        f"{TAXON_NAMES_CS_FILE}"
    )


# ============================================================
# 4. NAČÍST TAXONY Z AKTUÁLNÍHO PROJECTN
# ============================================================

if not COMPACT_NODES_FILE.exists():

    raise FileNotFoundError(
        "\nChybí soubor:\n"
        f"{COMPACT_NODES_FILE}\n\n"
        "Ten vytváří první průchod tree.py. "
        "Pipeline proto nemůže pokračovat."
    )


nodes_df = pd.read_csv(
    COMPACT_NODES_FILE,
    sep="\t"
)


if "taxid" not in nodes_df.columns:

    raise ValueError(
        f"{COMPACT_NODES_FILE} "
        "neobsahuje sloupec 'taxid'."
    )


wanted_taxids = sorted(
    set(
        nodes_df["taxid"]
        .dropna()
        .astype(int)
    )
)


print()
print(
    f"Taxonů v COMPACT stromu: "
    f"{len(wanted_taxids)}"
)


# ============================================================
# 5. NAČÍST OBECNOU TSV CACHE
# ============================================================

cache_df = pd.read_csv(
    TAXON_NAMES_CACHE_FILE,
    sep="\t",
    dtype={
        "taxid": "Int64",
        "lang": "string",
        "label": "string",
        "source": "string",
        "status": "string"
    }
)


required_cache_columns = {
    "taxid",
    "lang",
    "label",
    "source"
}


if not required_cache_columns.issubset(
    cache_df.columns
):

    raise ValueError(
        f"{TAXON_NAMES_CACHE_FILE} musí obsahovat alespoň sloupce: "
        "taxid, lang, label, source"
    )


# ------------------------------------------------------------
# ZPĚTNÁ KOMPATIBILITA SE STAROU CACHE
# ------------------------------------------------------------
#
# Starší cache neměla sloupec "status".
# U starých řádků tedy odvodíme:
#   neprázdný label -> found
#
# Nové řádky mohou mít:
#   status = found
#   status = not_found
#
# "not_found" znamená:
# Wikidata byla úspěšně dotázána, ale český label pro taxon nenašla.

if "status" not in cache_df.columns:

    cache_df["status"] = "found"


cache_df = cache_df.dropna(
    subset=[
        "taxid",
        "lang"
    ]
)


cache_df["lang"] = (
    cache_df["lang"]
    .astype(str)
    .str.strip()
)


# label smí být prázdný u status=not_found
cache_df["label"] = (
    cache_df["label"]
    .fillna("")
    .astype(str)
    .str.strip()
)


cache_df["status"] = (
    cache_df["status"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.lower()
)


# Staré řádky bez statusu / s prázdným statusem:
# pokud mají label, jsou "found".
cache_df.loc[
    (cache_df["status"] == "")
    & (cache_df["label"] != ""),
    "status"
] = "found"


cache_df = cache_df[
    cache_df["lang"] != ""
]


cache_cs_df = cache_df[
    cache_df["lang"] == LANGUAGE
]


# České názvy, které skutečně existují
cache_cs_found_df = cache_cs_df[
    (cache_cs_df["status"] == "found")
    & (cache_cs_df["label"] != "")
]


cached_cs_names = dict(
    zip(
        cache_cs_found_df["taxid"].astype(int),
        cache_cs_found_df["label"]
    )
)


# Taxony, na které už jsme se úspěšně ptali,
# ale Wikidata český label neměla.
cached_cs_not_found_taxids = set(
    cache_cs_df.loc[
        cache_cs_df["status"] == "not_found",
        "taxid"
    ]
    .astype(int)
)


print(
    f"Českých názvů v cache: "
    f"{len(cached_cs_names)}"
)

print(
    f"Taxonů v cache označených jako bez českého názvu: "
    f"{len(cached_cs_not_found_taxids)}"
)


# ============================================================
# 6. NAČÍST RUČNÍ ČESKÉ OPRAVY / OVERRIDE
# ============================================================

manual_cs_df = pd.read_csv(
    TAXON_NAMES_CS_FILE,
    sep="\t",
    dtype={
        "taxid": "Int64",
        "name_cs": "string"
    }
)


required_manual_columns = {
    "taxid",
    "name_cs"
}


if not required_manual_columns.issubset(
    manual_cs_df.columns
):

    raise ValueError(
        f"{TAXON_NAMES_CS_FILE} musí obsahovat "
        "sloupce 'taxid' a 'name_cs'."
    )


manual_cs_df = manual_cs_df.dropna(
    subset=[
        "taxid",
        "name_cs"
    ]
)


manual_cs_df["name_cs"] = (
    manual_cs_df["name_cs"]
    .astype(str)
    .str.strip()
)


manual_cs_df = manual_cs_df[
    manual_cs_df["name_cs"] != ""
]


manual_cs_names = dict(
    zip(
        manual_cs_df["taxid"].astype(int),
        manual_cs_df["name_cs"]
    )
)


print(
    f"Ručních / finálních českých názvů: "
    f"{len(manual_cs_names)}"
)


# ============================================================
# 7. URČIT, NA CO SE MUSÍME PTÁT WIKIDAT
# ============================================================

# "checked" znamená:
# - název už máme,
# - nebo už víme, že na Wikidatech český label není,
# - nebo je taxon ručně vyplněný.
#
# Díky tomu se na stejné "not_found" taxony neptáme při každém běhu.

checked_cs_taxids = (
    set(cached_cs_names)
    | set(cached_cs_not_found_taxids)
    | set(manual_cs_names)
)


missing_taxids = [
    taxid
    for taxid in wanted_taxids
    if taxid not in checked_cs_taxids
]


print(
    f"Taxonů, které je potřeba dohledat na Wikidatech: "
    f"{len(missing_taxids)}"
)


# ============================================================
# 8. WIKIDATA DOTAZ PRO JEDNU DÁVKU
# ============================================================

def fetch_labels_from_wikidata(
    taxids,
    language
):
    """
    Pro zadané NCBI Tax IDs najde na Wikidatech
    label v požadovaném jazyce.

    Při 429 nebo dočasné chybě několikrát zopakuje request.
    """

    values = " ".join(
        f'"{taxid}"'
        for taxid in taxids
    )


    query = f"""
    SELECT DISTINCT ?taxid ?label
    WHERE {{
        VALUES ?taxid {{
            {values}
        }}

        ?item wdt:P685 ?taxid .
        ?item rdfs:label ?label .

        FILTER(
            LANG(?label) = "{language}"
        )
    }}
    """


    last_error = None


    for attempt in range(MAX_RETRIES):

        try:

            response = requests.get(
                WIKIDATA_SPARQL_URL,
                params={
                    "query": query,
                    "format": "json"
                },
                headers=HEADERS,
                timeout=90
            )


            if response.status_code == 429:

                retry_after = (
                    response.headers.get(
                        "Retry-After"
                    )
                )

                if retry_after:
                    wait_seconds = float(
                        retry_after
                    )
                else:
                    wait_seconds = (
                        RETRY_BASE_SECONDS
                        * (2 ** attempt)
                    )

                print(
                    f"Wikidata hlásí 429. "
                    f"Čekám {wait_seconds:.0f} s "
                    f"a zkusím dávku znovu..."
                )

                time.sleep(
                    wait_seconds
                )

                continue


            response.raise_for_status()

            data = response.json()

            result = {}


            for binding in data[
                "results"
            ][
                "bindings"
            ]:

                taxid = int(
                    binding[
                        "taxid"
                    ][
                        "value"
                    ]
                )

                label = (
                    binding[
                        "label"
                    ][
                        "value"
                    ]
                    .strip()
                )


                if (
                    taxid not in result
                    and label != ""
                ):

                    result[taxid] = label


            return result


        except requests.RequestException as exc:

            last_error = exc

            wait_seconds = (
                RETRY_BASE_SECONDS
                * (2 ** attempt)
            )

            print(
                f"Wikidata request selhal: {exc}"
            )

            if attempt < MAX_RETRIES - 1:

                print(
                    f"Čekám {wait_seconds:.0f} s "
                    "a zkusím dávku znovu..."
                )

                time.sleep(
                    wait_seconds
                )


    if last_error is not None:
        raise last_error

    raise RuntimeError(
        "Wikidata request selhal i po opakovaných pokusech."
    )


# ============================================================
# 8B. ULOŽIT CHECKPOINT DO CACHE
# ============================================================

def save_batch_checkpoint(
    found_names,
    not_found_taxids
):
    """
    Okamžitě uloží výsledek jedné úspěšné Wikidata dávky do cache.

    Díky tomu se postup neztratí při:
    - Ctrl+C
    - pádu skriptu
    - výpadku internetu
    - vypnutí počítače

    Ukládá:
    - status=found pro nalezené české labely
    - status=not_found pro úspěšně zkontrolované taxony bez českého labelu
    """

    global cache_df

    rows = []


    for taxid, label in sorted(
        found_names.items()
    ):

        rows.append(
            {
                "taxid": taxid,
                "lang": LANGUAGE,
                "label": label,
                "source": "wikidata",
                "status": "found"
            }
        )


    for taxid in sorted(
        not_found_taxids
    ):

        # Pro jistotu: pokud taxon v téže dávce label měl,
        # nesmí se uložit jako not_found.
        if taxid in found_names:
            continue

        rows.append(
            {
                "taxid": taxid,
                "lang": LANGUAGE,
                "label": "",
                "source": "wikidata",
                "status": "not_found"
            }
        )


    if not rows:
        return


    batch_cache_df = pd.DataFrame(
        rows
    )


    cache_df = pd.concat(
        [
            cache_df,
            batch_cache_df
        ],
        ignore_index=True
    )


    # Jeden taxid + jazyk = jeden aktuální stav.
    # keep="last" dovolí později přepsat staré not_found nově nalezeným labelem.
    cache_df = cache_df.drop_duplicates(
        subset=[
            "taxid",
            "lang"
        ],
        keep="last"
    )


    cache_df = cache_df.sort_values(
        by=[
            "taxid",
            "lang"
        ]
    )


    # Zápis přes dočasný soubor:
    # nejdřív vytvoříme nový kompletní TSV a teprve potom nahradíme původní.
    # Je to bezpečnější při přerušení během zápisu.
    temp_cache_file = (
        TAXON_NAMES_CACHE_FILE
        .with_suffix(
            TAXON_NAMES_CACHE_FILE.suffix + ".tmp"
        )
    )


    cache_df.to_csv(
        temp_cache_file,
        sep="\t",
        index=False
    )


    temp_cache_file.replace(
        TAXON_NAMES_CACHE_FILE
    )

# ============================================================
# 9. DOTAZOVAT WIKIDATA POUZE NA CHYBĚJÍCÍ TAXONY
# ============================================================

new_names = {}
new_not_found_taxids = set()


if len(missing_taxids) == 0:

    print()
    print(
        "Cache pokrývá všechny potřebné taxony. "
        "Na Wikidata se neposílá žádný request."
    )

else:

    total_batches = (
        len(missing_taxids)
        + BATCH_SIZE
        - 1
    ) // BATCH_SIZE


    for batch_number, start in enumerate(
        range(
            0,
            len(missing_taxids),
            BATCH_SIZE
        ),
        start=1
    ):

        batch = missing_taxids[
            start:start + BATCH_SIZE
        ]


        print()
        print(
            f"Wikidata dávka "
            f"{batch_number}/{total_batches} "
            f"({len(batch)} taxonů)..."
        )


        batch_succeeded = False

        try:

            batch_names = (
                fetch_labels_from_wikidata(
                    batch,
                    LANGUAGE
                )
            )

            batch_succeeded = True

        except requests.RequestException as exc:

            print(
                "POZOR: Wikidata dotaz selhal "
                "i po opakovaných pokusech:"
            )

            print(exc)

            print(
                "Tuto dávku přeskakuji. "
                "Při příštím spuštění se zkusí znovu."
            )

            batch_names = {}


        print(
            f"Nalezeno českých názvů: "
            f"{len(batch_names)}"
        )


        # Nalezené labely
        for taxid, label in batch_names.items():

            if (
                taxid not in checked_cs_taxids
                and taxid not in new_names
            ):

                new_names[taxid] = label


        # DŮLEŽITÉ:
        # Jen pokud celý request proběhl úspěšně,
        # označíme taxony bez výsledku jako "not_found".
        #
        # Kdyby request selhal (429, timeout...), nesmíme z toho
        # usoudit, že český název neexistuje.
        if batch_succeeded:

            found_in_batch = set(
                batch_names
            )

            not_found_in_batch = (
                set(batch)
                - found_in_batch
            )

            new_not_found_taxids.update(
                not_found_in_batch
            )

            print(
                f"Bez českého labelu na Wikidatech: "
                f"{len(not_found_in_batch)}"
            )

            # CHECKPOINT:
            # úspěšnou dávku zapíšeme do cache HNED,
            # ne až na konci celého skriptu.
            save_batch_checkpoint(
                batch_names,
                not_found_in_batch
            )

            print(
                "Checkpoint této dávky uložen do cache."
            )


        if batch_number < total_batches:

            time.sleep(
                PAUSE_SECONDS
            )


# ============================================================
# 10. SOUHRN CHECKPOINTŮ
# ============================================================

# Všechny úspěšné dávky už byly zapsány do cache průběžně.
# Tady tedy nic dalšího nepřidáváme, jen vypíšeme souhrn.

print()
print(
    f"Nově nalezených českých názvů v tomto běhu: "
    f"{len(new_names)}"
)

print(
    f"Nově zapamatovaných taxonů bez českého názvu v tomto běhu: "
    f"{len(new_not_found_taxids)}"
)


# ============================================================
# 11. ZNOVU SESTAVIT ČESKÝ SOUBOR PRO TREE.PY
# ============================================================

# Priorita:
# 1. ruční hodnota v taxon_names_cs.tsv
# 2. cache z Wikidat

# cache_df už obsahuje i checkpointy z právě dokončených dávek.
current_cache_cs_found_df = cache_df[
    (cache_df["lang"] == LANGUAGE)
    & (cache_df["status"] == "found")
    & (cache_df["label"] != "")
]


final_cs_names = dict(
    zip(
        current_cache_cs_found_df["taxid"].astype(int),
        current_cache_cs_found_df["label"]
    )
)


for taxid, label in manual_cs_names.items():

    final_cs_names[
        taxid
    ] = label


output_cs_df = pd.DataFrame(
    [
        {
            "taxid": taxid,
            "name_cs": label
        }
        for taxid, label in sorted(
            final_cs_names.items()
        )
    ]
)


output_cs_df.to_csv(
    TAXON_NAMES_CS_FILE,
    sep="\t",
    index=False
)


# ============================================================
# 12. SOUHRN
# ============================================================

still_missing = [
    taxid
    for taxid in wanted_taxids
    if taxid not in final_cs_names
]


print()
print(
    "========================================"
)

print(
    f"Projekt: "
    f"{PROJECT_ID}"
)

print(
    f"Taxonů v aktuálním stromu: "
    f"{len(wanted_taxids)}"
)

print(
    f"Českých názvů dostupných po doběhnutí: "
    f"{len(final_cs_names)}"
)

print(
    f"Taxonů v aktuálním stromu bez českého názvu: "
    f"{len(still_missing)}"
)

print(
    f"Cache: "
    f"{TAXON_NAMES_CACHE_FILE}"
)

print(
    f"Český soubor pro tree.py: "
    f"{TAXON_NAMES_CS_FILE}"
)

print(
    "========================================"
)

print()
print(
    "Při dalším běhu se Wikidata budou dotazovat jen na Tax ID, "
    "která ještě nebyla zkontrolována. Každá úspěšná dávka se "
    "ukládá hned, takže postup přežije i Ctrl+C."
)
