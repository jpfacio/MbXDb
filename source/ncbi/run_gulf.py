from pathlib import Path
import argparse
import json
import subprocess
import resource
from time import perf_counter
from datetime import timedelta
import functions as f

# Control keys

build_jsonl = True
download_srrs = True
assembly = True
mapping = True
binning = True
metadata = True

# Source JSON (DOI -> SRR list) supplied by the --ncbi flag in mxbdb-build.sh

parser = argparse.ArgumentParser(
    description="Run the Gulf NCBI acquisition/assembly pipeline."
)
parser.add_argument(
    "json",
    nargs="?",
    type=Path,
    default=Path("support_files/mining/ncbi/gulf.json"),
    help="source JSON mapping publication DOIs to SRR runs",
)
parser.add_argument(
    "--srrs",
    help="comma-separated SRRs to process (default: every SRR in the JSON)",
)
parser.add_argument(
    "--data-dir",
    type=Path,
    default=Path("Data"),
    help="database root where bins and metadata.csv are written",
)
parser.add_argument(
    "--workdir",
    type=Path,
    default=Path("tmp/gulf_tmps"),
    help="working directory for downloads, assemblies and bins",
)
parser.add_argument(
    "--log-dir",
    type=Path,
    default=Path("log"),
    help="directory where pipeline logs are appended",
)
args = parser.parse_args()

source_json = args.json
data_dir = args.data_dir
log = args.log_dir

# Defining directories and files

out_gulf_dir = args.workdir
out_srr_dir = out_gulf_dir / "srr"
out_assembly_dir = out_gulf_dir / "assembly"
out_jsonl = out_srr_dir / f"{source_json.stem}.jsonl"
out_bin_dir = out_gulf_dir / "bins"
log_run = log / "run.log"

log.mkdir(parents=True, exist_ok=True)
out_gulf_dir.mkdir(parents=True, exist_ok=True)
out_srr_dir.mkdir(parents=True, exist_ok=True)
out_assembly_dir.mkdir(parents=True, exist_ok=True)
out_bin_dir.mkdir(parents=True, exist_ok=True)

with open(source_json) as fh:
    source = json.load(fh)

doi_by_srr = {
    srr: doi
    for doi, runs in source.items()
    for srr in runs
}

if args.srrs:
    requested = [s.strip() for s in args.srrs.split(",") if s.strip()]
    srrs = [s for s in requested if s in doi_by_srr]
    missing = [s for s in requested if s not in doi_by_srr]
    if missing:
        print(f"[WARN] SRRs not in {source_json.name}: {', '.join(missing)}")
else:
    srrs = list(doi_by_srr)

# Build gulf.jsonl from SRR metadata

# get_srr_info_batch does not carry DOIs, so they are stamped here from the
# source JSON at write time; keeps the metadata anchor's study_id populated.

if build_jsonl:

    fetch_start = perf_counter()

    srr_info = f.acq.get_srr_info_batch(srrs)

    for row in srr_info:
        row["doi"] = doi_by_srr.get(row.get("srr"), "")

    print(f"SRRs : {len(srr_info)}")

    with open(out_jsonl, "w") as fh:
        for row in srr_info:
            fh.write(f"{row}\n")

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", str(out_gulf_dir)],
        capture_output=True,
        text=True,
        check=True
    )

    project_size = space.stdout.split()[0]

    with open(log_run, 'a') as run:
        run.write(
            "#####  BUILD SRR METADATA CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Peak memory (children): {child_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n"
            f"SRR records: {len(srr_info)}\n\n"
        )

# Download reads for the recorded SRRs

if download_srrs:

    print('Downloading reads with fasterq-dump')

    fetch_start = perf_counter()

    completed = f.acq.download_srrs(
        srrs=srrs,
        outdir=str(out_srr_dir),
    )

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", str(out_srr_dir)],
        capture_output=True,
        text=True,
        check=True
    )

    project_size = space.stdout.split()[0]

    with open(log_run, 'a') as run:
        run.write(
            "#####  DOWNLOAD SRRS CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Peak memory (children): {child_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n"
            f"SRR completed: {completed:,}\n\n"
        )

# Assemble the reads of each run

if assembly:
    
    fetch_start = perf_counter()
    
    completed = f.ass.assembly(
        srrs=srrs,
        read_dir=out_srr_dir,
        outdir=out_assembly_dir
    )

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    
    space = subprocess.run(
        ["du", "-sh", str(out_assembly_dir)],
        capture_output=True,
        text=True,
        check=True
    )
    
    project_size = space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
        run.write(
            "#####  ASSEMBLE SRRS CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Peak memory (children): {child_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n"
            f"SRRs assembled: {completed:,}\n\n"
        )

if mapping:

    fetch_start = perf_counter()

    completed_idx = f.ass.build_index(out_assembly_dir)
    completed = f.ass.map_reads(
        srrs=srrs,
        read_dir=out_srr_dir,
        outdir=out_assembly_dir,
    )

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", str(out_assembly_dir)],
        capture_output=True,
        text=True,
        check=True
    )

    project_size = space.stdout.split()[0]

    with open(log_run, 'a') as run:
        run.write(
            "#####  MAP READS CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Peak memory (children): {child_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n"
            f"Indexes built: {completed_idx:,}\n"
            f"SRRs mapped: {completed:,}\n\n"
        )

if binning:

    fetch_start = perf_counter()

    completed = f.ass.bin_contigs(
        srrs=srrs,
        outdir=out_assembly_dir,
        bindir=out_bin_dir,
    )

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", str(out_bin_dir)],
        capture_output=True,
        text=True,
        check=True
    )

    project_size = space.stdout.split()[0]

    with open(log_run, 'a') as run:
        run.write(
            "#####  BIN CONTIGS CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Peak memory (children): {child_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n"
            f"SRRs binned: {completed:,}\n\n"
        )

if metadata:

    gulf_bins = data_dir / "Raw/Bins"
    out_metadata = data_dir / "metadata.csv"
    f.metadata.move_bins(out_bin_dir, gulf_bins)
    f.metadata.fetch_metadata(out_jsonl, gulf_bins, out_metadata)
