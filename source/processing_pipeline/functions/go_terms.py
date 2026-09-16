
import pandas as pd
import requests
import json
import re
import time
import ast
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from goatools.obo_parser import GODag

def create_database_metadata(genes: pd.DataFrame) -> pd.DataFrame:
    
    """This function defines a core element of the GO module, is extracts the dbs
    from the Genes entity and reorganizes the data according to databases origin.
    Then the function stores the ids in respective columns, returning a dataframe
    to be processed later in the pipeline. 
    
    The reason for this is the data format delivered by Bakta, which mix all the ids
    retrieved from databases in an unique string, the data needs to be saved in lists 
    to be processed easily.
    
    Args:
        genes (pd.DataFrame): Genes entity saved in the Entities module

    Returns:
        pd.DataFrame: Database metadata dataframe
    """
    databases = genes[['Gene_Tag', 'Databases']]
    databases = databases.dropna(subset=['Databases'])
    databases['Databases'] = databases['Databases'].str.split(', ')
    databases['Databases'] = databases['Databases'].apply(
        lambda x: [item for item in x if 'SO' not in item])
    
    for idx, row in databases.iterrows():
        for item in row["Databases"]:
            pref, id = map(str.strip, item.split(":", 1))

            if pref not in databases.columns:
                databases[pref] = None

            current = databases.loc[idx, pref]

            if pd.isna(current):
                databases.loc[idx, pref] = id
            else:
                databases.loc[idx, pref] = f"{current},{id}"
            
    if 'UniRef' in databases.columns:
        databases['UniRef'] = databases['UniRef'].str.split('_', expand=True)[1]
    
    return databases

def uniref2annotations(uniref: str, session: requests.Session, retries: int = 3):

    """This function queries the UniParc search endpoint once per UniRef
    representative accession and extracts both the UniParc IDs and the
    integrated InterPro entries from the same response. Every search hit
    carries its InterProScan matches in the sequenceFeatures field, where
    each feature holds the integrated InterPro group (e.g. IPR000001) when
    the member database signature has been integrated into an InterPro
    entry. Unintegrated signatures are skipped and the InterPro IDs of all
    hits are merged preserving first-seen order.

    Args:
        uniref (str): UniRef representative accession, e.g. A0A3M2J779
        session (requests.Session): Session to reuse the connection pool
        retries (int): Number of attempts before giving up

    Returns:
        tuple: (uniref, UniParc ID list or None, InterPro ID list or None)

    Raises:
        requests.exceptions.RequestException: after exhausting all retries,
        so callers can tell a real failure apart from a legitimate None
        ("no matches") and avoid caching failed lookups.
    """

    url = f"https://rest.uniprot.org/uniparc/search?query={uniref}"

    for attempt in range(retries):

        try:
            r = session.get(url, timeout=30)
            r.raise_for_status()

            data = r.json()

            results = data.get("results", [])

            uniparc = [
                x["uniParcId"]
                for x in results
            ]

            interpro = []
            for result in results:
                for feature in result.get("sequenceFeatures", []):
                    group = feature.get("interproGroup") or {}
                    if group.get("id"):
                        interpro.append(group["id"])

            interpro = list(dict.fromkeys(interpro))

            return uniref, uniparc or None, interpro or None

        except requests.exceptions.RequestException as e:

            if attempt < retries - 1:
                time.sleep(5)
            else:
                print(f"[ERROR] {uniref}: {e}")
                raise


def _load_uniref_cache(cache_dir) -> dict:

    """Loads previously resolved UniRef records from the cache directory.
    Each cache file is a JSON record {"uniref": ..., "uniparc": [...] or
    null, "interpro": [...] or null}.

    Args:
        cache_dir (str): Directory holding the per-ID result cache

    Returns:
        dict: {accession: {"uniparc": list or None,
                           "interpro": list or None}}
    """
    cache_dir = Path(cache_dir)
    cached = {}
    if cache_dir.is_dir():
        for f in cache_dir.glob("*.json"):
            with open(f) as fh:
                rec = json.load(fh)
            if "uniref" in rec:
                cached[rec["uniref"]] = {
                    "uniparc": rec.get("uniparc"),
                    "interpro": rec.get("interpro")
                }
    return cached


def _save_uniref_cache(cache_dir, uniref, uniparc, interpro) -> None:

    """Persists a single UniRef lookup result as a JSON file in the cache
    directory.

    Args:
        cache_dir (str): Directory holding the per-ID result cache
        uniref (str): UniRef representative accession
        uniparc (list): UniParc IDs found for the accession, or None
        interpro (list): InterPro IDs found for the accession, or None
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(cache_dir / f"{uniref}.json", "w") as fh:
        json.dump(
            {"uniref": uniref, "uniparc": uniparc, "interpro": interpro},
            fh
        )

def _parse_id_cell(ids) -> list:

    """Turns a column value holding identifier lists into a plain list of ID
    strings. Handles Python lists, the list reprs written to CSV by to_csv
    (e.g. "['UPI00003796B2']") and bare IDs.

    Args:
        ids: A single cell value from a list-valued column.

    Returns:
        list: ID strings extracted from the value.
    """
    if isinstance(ids, list):
        return [i for i in ids if isinstance(i, str)]
    if isinstance(ids, str):
        ids = ids.strip()
        if not ids:
            return []
        if ids.startswith("["):
            try:
                parsed = ast.literal_eval(ids)
            except (ValueError, SyntaxError):
                return []
            if isinstance(parsed, list):
                return [i for i in parsed if isinstance(i, str)]
            return []
        return [ids]
    return []

def fetch_uniref_annotations(data_meta: pd.DataFrame, max_workers: int = 20,
                             cache_dir: str = "tmp/uniref2ip_cache"):

    """This function maps every distinct UniRef representative accession in
    the database metadata to its UniParc entries and integrated InterPro
    entries, storing the results in new UniParc and InterPro columns. Both
    come from a single HTTP call per accession: the UniParc search response
    already carries the InterProScan matches of every hit, so no second
    round of per-entry calls is needed. Only integrated InterPro entries
    (IPR...) are kept, as a list of IDs per row.

    Results are cached per accession in `cache_dir` (one JSON file per ID),
    so an interrupted run can resume: already-resolved IDs are skipped and
    every completed lookup is written to disk as it finishes. Failed lookups
    are NOT cached, so they are retried on the next run.

    Args:
        data_meta (pd.DataFrame): Database metadata with a UniRef column
        max_workers (int): Number of parallel workers
        cache_dir (str): Directory holding the per-ID result cache

    Returns:
        pd.DataFrame: Database metadata with added UniParc/InterPro columns
    """

    query = sorted(
        data_meta["UniRef"]
        .dropna()
        .unique()
        .tolist()
    )

    cached = _load_uniref_cache(cache_dir)
    to_do = [acc for acc in query if acc not in cached]

    print(f"Processing {len(to_do):,} unique UniRef IDs "
          f"({len(cached):,} already cached)...\n")

    session = requests.Session()

    mapping = dict(cached)

    def worker(uniref):
        return uniref2annotations(uniref, session)

    if to_do:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:

            futures = {
                executor.submit(worker, acc): acc
                for acc in to_do
            }

            for future in tqdm(
                    as_completed(futures),
                    total=len(to_do),
                    desc="Mapping UniRef → UniParc/InterPro",
                    unit="IDs"
            ):
                try:
                    uniref, uniparc, interpro = future.result()
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")
                    continue
                mapping[uniref] = {"uniparc": uniparc, "interpro": interpro}
                _save_uniref_cache(cache_dir, uniref, uniparc, interpro)

    session.close()

    uniparc_map = {
        acc: rec["uniparc"] for acc, rec in mapping.items()
    }
    interpro_map = {
        acc: rec["interpro"] for acc, rec in mapping.items()
    }

    data_meta["UniParc"] = data_meta["UniRef"].map(uniparc_map)
    data_meta["InterPro"] = data_meta["UniRef"].map(interpro_map)

    mapped_upi = data_meta["UniParc"].notna().sum()
    mapped_ipr = data_meta["InterPro"].notna().sum()

    print("\nFinished!")
    print(f"Mapped UniParc {mapped_upi:,} / {len(data_meta):,} rows.")
    print(f"Mapped InterPro {mapped_ipr:,} / {len(data_meta):,} rows.")

    return data_meta


def parse_interpro2go(path: str) -> dict:

    """Parses the InterPro2GO mapping file into a dictionary mapping each
    InterPro entry (IPRxxxxxx) to its set of GO accessions.

    The file is available at:
        https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go
    and is line-oriented: comment lines start with "!" and data lines look like

        InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677

    Args:
        path (str): Path to the interpro2go file

    Returns:
        dict: {InterPro ID: set of GO accessions, e.g. {"GO:0003677"}}
    """
    ipr2go = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("!"):
                continue
            left, _, right = line.partition("> ")
            ipr = left.split(":", 1)[1].split(" ", 1)[0]
            go_ids = set(re.findall(r"GO:\d+", right))
            ipr2go.setdefault(ipr, set()).update(go_ids)
            
    return ipr2go


def fetch_go_terms(data_meta: pd.DataFrame, ipr2go: dict) -> pd.DataFrame:

    """Adds a GO column to the database metadata by mapping each gene's
    InterPro entries to GO terms through the InterPro2GO mapping.

    Args:
        data_meta (pd.DataFrame): Database metadata with an InterPro column
        ipr2go (dict): InterPro -> set of GO accessions (see parse_interpro2go)

    Returns:
        pd.DataFrame: Database metadata with the added GO column
    """
    def go_for(iprs):
        go = set()
        for ipr in _parse_id_cell(iprs):
            go.update(ipr2go.get(ipr, set()))
        return sorted(go) or None

    data_meta["GO"] = data_meta["InterPro"].apply(go_for)

    mapped = data_meta["GO"].notna().sum()

    print(f"Mapped {mapped:,} / {len(data_meta):,} rows.")

    return data_meta


def resolve_go_names(data_meta: pd.DataFrame, obo_path: str) -> pd.DataFrame:

    """Adds a GO_Name column by mapping each GO accession in the GO column
    to its human-readable name using the GO OBO ontology.

    Args:
        data_meta (pd.DataFrame): Database metadata with a GO column
        obo_path (str): Path to the GO OBO file (e.g. go-basic.obo)

    Returns:
        pd.DataFrame: Database metadata with the added GO_Name column
    """
    godag = GODag(obo_path)

    def names_for(go_ids):
        ids = _parse_id_cell(go_ids)
        if not ids:
            return None
        names = sorted(godag[go_id].name for go_id in ids if go_id in godag)
        return names or None

    data_meta["GO_Name"] = data_meta["GO"].apply(names_for)

    mapped = data_meta["GO_Name"].notna().sum()
    print(f"Resolved {mapped:,} / {len(data_meta):,} rows to GO names.")

    return data_meta

