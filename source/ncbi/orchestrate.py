"""Plan and run NCBI acquisitions in storage-friendly chunks.

The source JSONs are read in the order given by --sources, preserving the DOI
and SRR order inside each file. Runs already present in Data/metadata.csv are
skipped. Depending on the mode:

  chunks N    process the remaining runs N at a time until none are left
  isolated V  process one batch from a {doi: [srr]} JSON object or JSON file

Every batch runs the matching source script and the processing pipeline inside
a staging Data tree. The staging tables and raw files are merged into the
master Data tree and the batch tmp directory is removed afterwards. Each batch
also gets its own log directory under log/chunks/.
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run NCBI sources in chunks or as an isolated batch."
    )
    parser.add_argument("--mode", choices=["chunks", "isolated"], required=True)
    parser.add_argument("--value", required=True,
                        help="chunk size or isolated-run JSON value/path")
    parser.add_argument("--sources", required=True,
                        help="comma-separated source names, in order")
    parser.add_argument("--json-dir", type=Path,
                        default=Path("support_files/mining/ncbi"))
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--tmp-root", type=Path, default=Path("tmp"))
    parser.add_argument("--log-dir", type=Path, default=Path("log"))
    parser.add_argument("--bakta-db", type=Path, default=None)
    parser.add_argument("--keep-tmps", action="store_true")
    return parser.parse_args()


def processed_srrs(metadata: Path) -> set[str]:
    if not metadata.exists():
        return set()

    with open(metadata, newline="") as fh:
        reader = csv.DictReader(fh)
        return {row["srr"] for row in reader if row.get("srr")}


def load_ordered(sources: list[str], json_dir: Path) -> list[tuple[str, str, Path]]:
    ordered = []
    seen = set()

    for source in sources:
        source_json = json_dir / f"{source}.json"
        with open(source_json) as fh:
            data = json.load(fh)

        for runs in data.values():
            for srr in runs:
                srr = srr.strip()
                if not srr or srr in seen:
                    continue
                seen.add(srr)
                ordered.append((srr, source, source_json))

    return ordered


def load_isolated(value: str, source: str,
                  tmp_root: Path) -> list[tuple[str, str, Path]]:
    path = Path(value)

    if path.is_file():
        try:
            with open(path) as fh:
                data = json.load(fh)
        except json.JSONDecodeError as exc:
            raise SystemExit(f"error: invalid JSON in {path}: {exc}")
    else:
        try:
            data = json.loads(value)
        except json.JSONDecodeError:
            raise SystemExit(
                "error: --isolated-run expects a {doi: [srr]} JSON object "
                "or a path to a JSON file"
            )

    if not isinstance(data, dict) or not data:
        raise SystemExit("error: --isolated-run must be a non-empty JSON object")

    normalized = {}
    ordered = []
    seen = set()

    for raw_doi, runs in data.items():
        doi = str(raw_doi).strip()
        if not doi:
            raise SystemExit("error: --isolated-run DOI keys must be non-empty")
        if not isinstance(runs, list):
            raise SystemExit(f"error: DOI '{doi}' must map to a list of SRRs")

        cleaned = []
        for raw_srr in runs:
            srr = str(raw_srr).strip()
            if not srr:
                continue
            cleaned.append(srr)
            if srr not in seen:
                seen.add(srr)
                ordered.append(srr)

        if cleaned:
            normalized[doi] = cleaned

    if not ordered:
        raise SystemExit("error: --isolated-run contains no SRRs")

    source_json = tmp_root / "isolated_input.json"
    source_json.parent.mkdir(parents=True, exist_ok=True)
    with open(source_json, "w") as fh:
        json.dump(normalized, fh, indent=2)

    return [(srr, source, source_json) for srr in ordered]


def positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise SystemExit(f"error: expected a positive integer, got '{value}'")
    if number <= 0:
        raise SystemExit("error: expected a positive integer")
    return number


def task_chunks(runs, size: int):
    for start in range(0, len(runs), size):
        yield runs[start:start + size]


def run_batch(chunk, index, args) -> None:
    chunk_root = args.tmp_root / f"ncbi_chunk_{index:04d}_{chunk[0][0]}"
    chunk_log = args.log_dir / "chunks" / chunk_root.name
    staging = chunk_root / "Data"
    workdir = chunk_root / "work"

    if chunk_root.exists():
        shutil.rmtree(chunk_root)

    for directory in (
        staging / "Raw/Bins",
        staging / "Raw/Processed",
        staging / "Entities",
        staging / "Reports",
        workdir,
        chunk_log,
        args.log_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    srrs = [srr for srr, _, _ in chunk]
    sources = sorted({source for _, source, _ in chunk})

    with open(args.log_dir / "run.log", "a") as index_log:
        index_log.write(
            f"Chunk {index}: {len(chunk)} SRR(s) -> {chunk_log}\n"
        )

    with open(chunk_log / "run.log", "a") as chunk_run:
        chunk_run.write(
            f"##### CHUNK {chunk_root.name} #####\n"
            f"Source(s): {', '.join(sources)}\n"
            f"SRRs: {', '.join(srrs)}\n\n"
        )

    grouped = OrderedDict()
    for srr, source, source_json in chunk:
        group = grouped.setdefault(source, {"json": source_json, "srrs": []})
        group["srrs"].append(srr)

    for source, group in grouped.items():
        run_script = Path("source/ncbi") / f"run_{source}.py"
        print(f"[{index}] {source}: {', '.join(group['srrs'])}")
        subprocess.run(
            [
                sys.executable, str(run_script), str(group["json"]),
                "--srrs", ",".join(group["srrs"]),
                "--data-dir", str(staging),
                "--workdir", str(workdir),
                "--log-dir", str(chunk_log),
            ],
            check=True,
        )

    pp_cmd = [
        sys.executable, "source/processing_pipeline/run_pp.py",
        "--data-dir", str(staging),
        "--tmp-dir", str(workdir),
        "--log-dir", str(chunk_log),
    ]
    if args.bakta_db is not None:
        pp_cmd += ["--bakta-db", str(args.bakta_db)]

    subprocess.run(pp_cmd, check=True)

    subprocess.run(
        [
            sys.executable, "source/processing_pipeline/merge_tables.py",
            "--src", str(staging),
            "--dest", str(args.data_dir),
        ],
        check=True,
    )

    if args.keep_tmps:
        print(f"[{index}] keeping {chunk_root}")
    else:
        shutil.rmtree(chunk_root)
        print(f"[{index}] removed {chunk_root}")


def main() -> None:
    args = parse_args()
    sources = [s.strip() for s in args.sources.split(",") if s.strip()]

    if not sources:
        raise SystemExit("error: no sources selected")

    done = processed_srrs(args.data_dir / "metadata.csv")

    if args.mode == "chunks":
        size = positive_int(args.value)
        ordered = load_ordered(sources, args.json_dir)
        runs = [item for item in ordered if item[0] not in done]
        batches = task_chunks(runs, size)
    else:
        if len(sources) != 1:
            raise SystemExit(
                "error: --isolated-run requires exactly one source"
            )
        ordered = load_isolated(args.value, sources[0], args.tmp_root)
        runs = [item for item in ordered if item[0] not in done]
        batches = [runs]

    if not runs:
        print("Nothing to do: every requested SRR is already in metadata.csv.")
        return

    for index, chunk in enumerate(batches, start=1):
        run_batch(chunk, index, args)

    print(f"Finished {len(runs):,} SRRs.")


if __name__ == "__main__":
    main()
