import ast
import csv
import gzip
import shutil
import subprocess
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from xml.etree import ElementTree

from tqdm import tqdm


def _natural_key(name: str) -> tuple:
    """Split a bin stem into (prefix, number) for numeric-aware ordering."""
    prefix, _, tail = name.partition(".")
    try:
        return prefix, int(tail)
    except ValueError:
        return name, 0


def _srr_of(stem: str) -> str:
    """SRR run prefix of a flat bin name (SRR21147275.1.fa -> SRR21147275)."""
    return stem.partition(".")[0]


def _bin_id(stem: str) -> str:
    """Bin identifier of a flat bin name, dropping the fasta suffix.

    Matches the name Bakta prints on its output ('SRR21147275.1.fa.gz' ->
    'SRR21147275.1'), so bins.csv.Bin joins against genes.csv.Bin.
    """
    return stem.removesuffix(".fa")


def _parse_biosample_attributes(xml_text: str) -> dict:
    """Extract attribute_name -> value pairs from eutils BioSample XML.

    Args:
        xml_text (str): raw XML returned by efetch db=biosample

    Returns:
        dict: BioSample attribute values (empty dict when the record has none)
    """
    root = ElementTree.fromstring(xml_text)
    attrs = {}
    for attrs_node in root.iter("Attributes"):
        for attribute in attrs_node.findall("Attribute"):
            name = attribute.get("attribute_name")
            if name:
                attrs[name] = (attribute.text or "").strip()
    return attrs


def _fetch_biosample(accession: str) -> dict:
    """Fetch BioSample attributes for a SAMN accession.

    Mirrors get_srr_info (acquisition.py): wget on NCBI eutils, a short
    sleep before every call to respect the ~3 requests/sec API cap, and up
    to five attempts with exponential backoff on failure or empty body.

    Args:
        accession (str): BioSample accession (e.g. "SAMN05791315")

    Returns:
        dict: attribute_name -> value for the accession

    Raises:
        RuntimeError: when every attempt fails or returns an empty body
    """
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=biosample&id={accession}&retmode=xml")

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
            return _parse_biosample_attributes(result.stdout)

        time.sleep(2 ** attempt)

    raise RuntimeError(f"biosample fetch failed for {accession}")


def _fetch_biosample_batch(accessions: list, max_workers: int = 4) -> dict:
    """Fetch BioSample attributes for many accessions in parallel.

    Failed lookups are reported and treated as empty attribute sets, so the
    metadata table still lists every bin; the gap stays visible as empty
    coord/date/depth fields rather than silently dropping samples.

    Args:
        accessions (list): SAMN accessions
        max_workers (int): number of concurrent fetchers

    Returns:
        dict: accession -> attribute dict
    """
    out = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_fetch_biosample, acc): acc
            for acc in accessions
        }

        for future in tqdm(
                as_completed(futures),
                total=len(accessions),
                desc="BioSample attributes",
                unit="samples"
        ):
            acc = futures[future]
            try:
                out[acc] = future.result()
            except Exception as e:
                out[acc] = {}
                print(f"[ERROR] {acc}: {e}")

    return out


def fetch_metadata(
    jsonl: Path = Path("tmp/gulf_tmps/srr/gulf.jsonl"),
    bindir: Path = Path("Data/Raw/Bins"),
    out_csv: Path = Path("tmp/metadata.csv"),
    max_workers: int = 4,
) -> int:
    """Build the metadata.csv anchor for the processing pipeline.

    All-in-one filler for the Gulf of Mexico subproject: every bin in the
    QC input directory is stored as its own row, the SRR base name goes into
    the srr column, and sample/project/study_id come straight from
    gulf.jsonl (biosample, bioproject, doi). The coord, date and depth
    columns are pulled live from the NCBI BioSample record of each SAMN
    accession via eutils (see _fetch_biosample).

    Bins are expected as flat *.fa.gz files (the shape move_bins produces
    and the QC stage consumes); the SRR of each bin is derived from its name
    (SRR21147275.1 -> SRR21147275). Bins not present in gulf.jsonl are
    skipped, keeping the table scoped to the gulf subproject.

    Args:
        jsonl (Path): gulf SRR metadata dump (one dict per line)
        bindir (Path): directory holding the flat gzipped bins
        out_csv (Path): where to write the metadata table
        max_workers (int): concurrent BioSample fetches

    Returns:
        int: number of metadata rows written
    """
    jsonl = Path(jsonl)
    bindir = Path(bindir)
    out_csv = Path(out_csv)

    by_srr = {}
    with open(jsonl) as fh:
        for line in fh:
            row = ast.literal_eval(line)
            if row.get("srr"):
                by_srr[row["srr"]] = row

    bins = sorted(
        bindir.glob("*.fa.gz"),
        key=lambda p: _natural_key(_bin_id(p.stem))
    )

    bin_srrs = sorted({_srr_of(_bin_id(p.stem)) for p in bins})

    accessions = sorted({
        by_srr[srr]["biosample"] for srr in bin_srrs
        if srr in by_srr and by_srr[srr].get("biosample")
    })
    attrs = _fetch_biosample_batch(accessions, max_workers=max_workers)

    columns = ["bin", "srr", "sample", "project", "study_id", "date", "coord", "depth"]
    rows = []

    for bin_path in bins:
        bin_id = _bin_id(bin_path.stem)
        srr = _srr_of(bin_id)
        if srr not in by_srr:
            continue

        row = by_srr[srr]
        biosample = row.get("biosample") or ""
        sample_attrs = attrs.get(biosample, {})

        rows.append([
            bin_id,
            srr,
            biosample,
            row.get("bioproject") or "",
            row.get("doi") or "",
            sample_attrs.get("collection_date") or "",
            sample_attrs.get("lat_lon") or "",
            sample_attrs.get("depth") or "",
        ])

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        writer.writerows(rows)

    print(f"Wrote {len(rows):,} bins "
          f"({len(set(r[1] for r in rows))} SRRs) to {out_csv}")
    return len(rows)


def move_bins(
    bindir: Path = Path("tmp/gulf_tmps/bins"),
    destdir: Path = Path("Data/Raw/Bins"),
) -> int:
    """Move gulf bins into the QC input directory, flattened and gzipped.

    The processing pipeline's QC stage expects a flat directory of *.fa.gz
    files (sum_stats.seqkit_summary globs Data/Raw/Bins/*.fa.gz). Gulf bins
    are produced per-SRR under tmp/gulf_tmps/bins, so each bin fasta is
    moved out of its SRR folder and gzipped while landing. Destinations that
    already exist are left untouched, making repeated runs harmless.

    Args:
        bindir (Path): gulf binning output (one subdirectory per SRR)
        destdir (Path): QC input directory (Data/Raw/Bins)

    Returns:
        int: number of bins moved
    """
    bindir = Path(bindir)
    destdir = Path(destdir)
    destdir.mkdir(parents=True, exist_ok=True)

    moved = 0
    for src_path in sorted(bindir.glob("*/[!.]*.fa")):
        dest_path = destdir / (src_path.name + ".gz")
        if dest_path.exists():
            continue

        with open(src_path, "rb") as fh_in, gzip.open(dest_path, "wb") as fh_out:
            shutil.copyfileobj(fh_in, fh_out)

        src_path.unlink()
        moved += 1

    print(f"Moved {moved:,} bins to {destdir} (gzipped)")
    return moved