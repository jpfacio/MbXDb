import argparse
from pathlib import Path
import subprocess
import resource
import json
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

# Defining directories and files

tmp = Path("tmp")
log = Path("log")
out_ncbi_dir = tmp / "ncbi_tmps"
out_srr_dir = out_ncbi_dir / "srr"
out_assembly_dir = out_ncbi_dir / "assembly"
out_jsonl = tmp / "ncbi.jsonl"
out_bin_dir = out_ncbi_dir / "bins"
log_run = log / "run.log"
json_dir = Path("support_files/ncbi")
ncbi_bins = Path("Data/Raw/Bins")
out_metadata = Path("Data/metadata.csv")

# Arguments

parser = argparse.ArgumentParser(
    prog="run_ncbi.py",
    description="Build the NCBI SRR stage from the jsons in support_files/ncbi."
)
parser.add_argument(
    "names",
    nargs="*",
    help="json names in support_files/ncbi (default: every json)"
)
parser.add_argument(
    "--chunks",
    type=int,
    default=None,
    help="process this many SRRs per chunk, clearing tmp between chunks"
)
args = parser.parse_args()

chunk_size = args.chunks

if chunk_size is not None and chunk_size < 1:
    raise SystemExit("--chunks must be a positive integer")

tmp.mkdir(parents=True, exist_ok=True)
log.mkdir(parents=True, exist_ok=True)
out_ncbi_dir.mkdir(parents=True, exist_ok=True)
out_srr_dir.mkdir(parents=True, exist_ok=True)
out_assembly_dir.mkdir(parents=True, exist_ok=True)
out_bin_dir.mkdir(parents=True, exist_ok=True)

# Resolve the NCBI inputs: explicit json names, or every json when none given

if args.names:
    json_files = [
        json_dir / (n if n.endswith(".json") else f"{n}.json")
        for n in args.names
    ]
else:
    json_files = sorted(json_dir.glob("*.json"))

missing = [str(p) for p in json_files if not p.exists()]
if missing:
    raise SystemExit(f"Unknown NCBI input(s): {', '.join(missing)}")

# Each json maps a publication DOI to the SRRs that belong to it. Collect the
# SRR pool and remember which DOI each run belongs to.

srrs = []
doi_by_srr = {}

for path in json_files:
    with open(path) as fh:
        project = json.load(fh)
    for doi, project_srrs in project.items():
        for srr in project_srrs:
            srrs.append(srr)
            doi_by_srr[srr] = doi

srrs = sorted(set(srrs))

print(f"NCBI inputs: {', '.join(p.stem for p in json_files)}")
print(f"SRRs: {len(srrs)}")

# Build ncbi.jsonl from SRR metadata

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
        ["du", "-sh", str(out_ncbi_dir)],
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

# Split the pool into chunks (a single chunk when --chunks is not given)

if chunk_size:
    chunks = [srrs[i:i + chunk_size] for i in range(0, len(srrs), chunk_size)]
else:
    chunks = [srrs]

print(f"Chunks: {len(chunks)}")

total_bins = 0

for index, chunk in enumerate(chunks, 1):

    if chunk_size:
        print(f"#####  CHUNK {index}/{len(chunks)} ({len(chunk)} SRRs)  #####")

        with open(log_run, 'a') as run:
            run.write(
                f"##########  CHUNK {index}/{len(chunks)} "
                f"({len(chunk)} SRRs)  ##########\n\n"
            )

    # Download reads for the chunk

    if download_srrs:

        print('Downloading reads with fasterq-dump')

        fetch_start = perf_counter()

        completed = f.acq.download_srrs(
            srrs=chunk,
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
            srrs=chunk,
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
            srrs=chunk,
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
            srrs=chunk,
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

    # Move the chunk's bins out of tmp before the working directories are cleared

    if metadata:
        total_bins += f.metadata.move_bins(out_bin_dir, ncbi_bins)

    # Free the per-chunk working directories, keeping ncbi.jsonl (in tmp/)

    if chunk_size and index < len(chunks):
        print("Clearing tmp for next chunk")
        f.util.clear_working(out_ncbi_dir)

# Regenerate the metadata anchor once, from every bin accumulated across chunks

if metadata:
    f.metadata.fetch_metadata(out_jsonl, ncbi_bins, out_metadata)

print(f"Bins produced: {total_bins:,}")

if metadata and total_bins == 0:
    print("[WARN] No bins were produced from the selected NCBI inputs.")
