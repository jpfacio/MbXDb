from pathlib import Path
import subprocess
import resource
import ast
from time import perf_counter
from datetime import timedelta
import functions as f

# Control keys

build_jsonl = True
download_srrs = True

# Defining directories and files

tmp = Path("tmp")
out_dir = tmp / "gulf_tmps"
out_dir.mkdir(exist_ok=True)
out = out_dir / "gulf.jsonl"
log = Path("log")
log_run = Path("log/run.log")

# Build gulf.jsonl from SRR metadata

if build_jsonl:

    srrs = ["SRR4342129", "SRR4342130", "SRR4342133", "SRR4342134", "SRR4342135", "SRR4342136", 
            "SRR21147274", "SRR21147275", "SRR21147276", "SRR21147277", "SRR21147278",
            "SRR21147279", "SRR21147280", "SRR21147281", "SRR21147282", "SRR21147283"]

    fetch_start = perf_counter()

    srr_info = f.acq.get_srr_info_batch(srrs)

    print(f"SRRs : {len(srr_info)}")

    with open(out, "w") as fh:
        for row in srr_info:
            fh.write(f"{row}\n")

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", '.'],
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
            f"Project size: {project_size}\n"
            f"SRR records: {len(srr_info)}\n\n"
        )

# ======================================================================
# TEST SUBSET ##########################################################
# Read gulf.jsonl and keep only ONE SRR run from the pool, to test the
# pipeline on a single run. The subset lives in memory only (test_srrs).
# REMOVE THIS WHOLE BLOCK (and the `srrs=test_srrs` argument in the
# download stage below) to download the full SRR pool from gulf.jsonl.
# ======================================================================

all_srrs = []
with open(out) as fh:
    for line in fh:
        all_srrs.append(ast.literal_eval(line))

test_srrs = [all_srrs[0]["srr"]]          # only the first SRR of the pool

# ======================================================================

# Download reads for the recorded SRRs

if download_srrs:

    print('Downloading reads with fasterq-dump')

    fetch_start = perf_counter()

    completed = f.acq.download_srrs(
        srrs=test_srrs,                   # TEST SUBSET: single SRR run only
        outdir=str(out_dir),
    )
    # Full pool instead:
    # completed = f.acq.download_srrs(jsonl=str(out), outdir=str(out_dir))

    fetch_elapsed = perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    space = subprocess.run(
        ["du", "-sh", str(out_dir)],
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
            f"Project size: {project_size}\n"
            f"SRR completed: {completed:,}\n\n"
        )