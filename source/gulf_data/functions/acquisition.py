import subprocess
import json
import csv
import time
import ast
from io import StringIO
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def get_gcf_info(accession: str) -> list:
    """Fetch the assembled MAG (GCA) summary for a BioProject.

    Used when the BioProject contains assembled metagenomes. The reads are
    already assembled, so the run field is left empty.

    Args:
        accession (str): BioProject accession (PRJNA...)

    Returns:
        list: One dict per MAG (accession, organism, bioproject, biosample, srr)
    """

    result = subprocess.run(
        ["datasets", "summary", "genome", "accession", accession,
         "--mag", "only", "--as-json-lines"],
        capture_output=True, text=True, check=True
    )

    rows = [json.loads(line) for line in result.stdout.splitlines()]

    subset = [
        {
            "accession": d["accession"],
            "organism": d["organism"]["organism_name"],
            "bioproject": d["assembly_info"]["bioproject_accession"],
            "biosample": d["assembly_info"]["biosample"]["accession"],
            "srr": None,
        }
        for d in rows
    ]

    return subset


def get_srr_info(accession: str) -> list:
    """Fetch eutils runinfo for an SRR run.

    Used when the BioProject only provides raw reads (no assembled metagenomes).
    The assembly accession is left empty until the reads get assembled.

    eutils caps ~3 requests/sec per IP without an API key; back-to-back calls
    trip it and wget exits 8 with an empty body. Each call paces itself and
    retries with backoff on failure or empty response.

    Args:
        accession (str): SRR accession (e.g. "SRR21147283")

    Returns:
        list: One dict per run (accession, organism, bioproject, biosample, srr)
    """

    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=sra&id={accession}&rettype=runinfo&retmode=text")

    for attempt in range(5):
        time.sleep(0.4)
        try:
            result = subprocess.run(
                ["wget", "-qO-", url],
                capture_output=True, text=True, check=True, timeout=30
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            time.sleep(2 ** attempt)
            continue

        if result.stdout.strip():
            return _parse_runinfo(result.stdout)

        time.sleep(2 ** attempt)

    raise RuntimeError(f"runinfo failed for {accession}")


def get_srr_info_batch(srrs: list, max_workers: int = 4) -> list:
    """Fetch runinfo for many SRRs in parallel.

    Mirrors the other modules' ThreadPoolExecutor + as_completed + tqdm
    pattern; pacing and retry live inside get_srr_info, so throttled calls
    self-heal. Failed lookups are reported and skipped.

    Args:
        srrs (list): SRR accessions
        max_workers (int): Number of concurrent fetchers

    Returns:
        list: flattened run rows from all SRRs
    """

    rows = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:

        futures = {
            executor.submit(get_srr_info, srr): srr
            for srr in srrs
        }

        for future in tqdm(
                as_completed(futures),
                total=len(srrs),
                desc="SRR runinfo",
                unit="runs"
        ):
            try:
                rows.extend(future.result())
            except Exception as e:
                print(f"[ERROR] {futures[future]}: {e}")

    return rows


def merge(gcf_rows: list, srr_rows: list) -> list:
    """Merge GCF and SRR rows into one list of dicts.

    Both inputs share the same structure (accession, organism, bioproject,
    biosample, srr). Every GCF assembly is kept as its own row; an SRR row is
    attached to a GCF row when it shares its biosample, and SRR rows with no
    matching assembly are emitted on their own.

    Args:
        gcf_rows (list): rows from get_gcf_info (srr empty)
        srr_rows (list): rows from get_srr_info (accession empty)

    Returns:
        list: one dict per GCF assembly, plus unmatched SRR rows
    """

    merged = {row["accession"]: dict(row) for row in gcf_rows}
    by_biosample = {row["biosample"]: row for row in merged.values()}

    for row in srr_rows:
        hit = by_biosample.get(row["biosample"])
        if hit is not None and hit["srr"] is None:
            hit["srr"] = row["srr"]
            if row["organism"]:
                hit["organism"] = row["organism"]
            continue

        key = row["accession"] or row["srr"]
        if key not in merged:
            merged[key] = dict(row)

    return list(merged.values())


def _parse_runinfo(text: str) -> list:
    """Parse eutils runinfo CSV text into the run-dict schema."""

    if not text.strip():
        return []

    rows = list(csv.DictReader(StringIO(text)))

    return [
        {
            "accession": None,
            "organism": r.get("ScientificName"),
            "bioproject": r.get("BioProject"),
            "biosample": r.get("BioSample"),
            "srr": r.get("Run"),
        }
        for r in rows
    ]


def download_srrs(
    jsonl: str = "tmp/gulf_tmps/gulf.jsonl",
    srrs: list = None,
    outdir: str = "tmp/gulf_tmps",
    max_workers: int = 3
) -> int:
    """Download FASTQ reads for a set of SRRs.

    Runs fasterq-dump in parallel (ThreadPoolExecutor + tqdm, mirroring the
    other modules). Runs already present in outdir are skipped, so restarted
    runs resume instead of re-downloading. Logging is left to the caller.

    Args:
        jsonl (str): path to the gulf output dump (one dict per line), read
            when srrs is not given
        srrs (list): explicit SRR accessions; when given, jsonl is not read
        outdir (str): directory where fasterq-dump writes the FASTQ files
        max_workers (int): number of concurrent fasterq-dump processes

    Returns:
        int: number of SRRs successfully downloaded
    """

    if srrs is not None:
        srrs = sorted(set(x for x in srrs if x))
    else:
        srrs = []
        with open(jsonl) as fh:
            for line in fh:
                row = ast.literal_eval(line)
                if row.get("srr"):
                    srrs.append(row["srr"])
        srrs = sorted(set(srrs))

    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    to_do = [s for s in srrs if not list(out_dir.glob(f"{s}*.fastq"))]

    print(f"Downloading {len(to_do):,} SRRs "
          f"({len(srrs) - len(to_do):,} already present)...\n")

    def worker(srr):
        subprocess.run(
            ["fasterq-dump", srr, "-O", str(out_dir), "-t", str(out_dir),
             "--split-files", "--progress", "-e", "4"],
            check=True
        )
        return srr

    completed = 0
    if to_do:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:

            futures = {
                executor.submit(worker, srr): srr
                for srr in to_do
            }

            for future in tqdm(
                    as_completed(futures),
                    total=len(to_do),
                    desc="fasterq-dump",
                    unit="runs"
            ):
                try:
                    future.result()
                    completed += 1
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")

    return completed