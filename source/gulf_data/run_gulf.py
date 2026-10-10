from pathlib import Path
import subprocess
import resource
import ast
import random
import shutil
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
out_gulf_dir = tmp / "gulf_tmps"
out_srr_dir = out_gulf_dir / "srr"
out_assembly_dir = out_gulf_dir / "assembly"
out_jsonl = out_srr_dir / "gulf.jsonl"
out_bin_dir = out_gulf_dir / "bins"
log_run = log / "run.log"

tmp.mkdir(parents=True, exist_ok=True)
log.mkdir(parents=True, exist_ok=True)
out_gulf_dir.mkdir(parents=True, exist_ok=True)
out_srr_dir.mkdir(parents=True, exist_ok=True)
out_assembly_dir.mkdir(parents=True, exist_ok=True)
out_bin_dir.mkdir(parents=True, exist_ok=True)

srrs = ["SRR4342129", "SRR4342130", "SRR4342133", "SRR4342134", "SRR4342135", "SRR4342136", 
            "SRR21147274", "SRR21147275", "SRR21147276", "SRR21147277", "SRR21147278",
            "SRR21147279", "SRR21147280", "SRR21147281", "SRR21147282", "SRR21147283"]

# Build gulf.jsonl from SRR metadata

# Hardcoded BioProject -> publication DOI. get_srr_info_batch does not carry
# DOIs, so they are stamped here at write time; keeps the metadata anchor's
# study_id populated across rebuilds.

bioproject_dois = {
    "PRJNA340003": "10.1038/s41597-025-05736-9",
    "PRJNA870083": "10.1128/aem.00799-26",
}

if build_jsonl:

    fetch_start = perf_counter()

    srr_info = f.acq.get_srr_info_batch(srrs)

    for row in srr_info:
        row["doi"] = bioproject_dois.get(row.get("bioproject"), "")

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

# ======================================================================
# TEST SUBSET ##########################################################
# Read gulf.jsonl and keep TWO random SRR runs from the pool, to test the
# pipeline on a small subset. The subset lives in memory only (test_srrs).
# REMOVE THIS WHOLE BLOCK (and the `srrs=test_srrs` argument in the
# download and assembly stages below) to run the full SRR pool from
# gulf.jsonl. Only read when a stage that uses it is enabled, otherwise
# a fresh checkout dies here before reaching any gated stage.
# ======================================================================

if download_srrs or assembly or mapping or binning:
    rng = random.Random(42)         

    all_srrs = []
    with open(out_jsonl) as fh:
        for line in fh:
            all_srrs.append(ast.literal_eval(line))

    pool = sorted({row["srr"] for row in all_srrs if row.get("srr")})

    k = min(2, len(pool))
    test_srrs = rng.sample(pool, k) if k else []

# ======================================================================

# Download reads for the recorded SRRs

if download_srrs:

    print('Downloading reads with fasterq-dump')

    fetch_start = perf_counter()

    completed = f.acq.download_srrs(
        srrs=test_srrs,                   # TEST SUBSET: two random SRR runs only
        outdir=str(out_srr_dir),
    )
    # Full pool instead:
    # completed = f.acq.download_srrs(jsonl=str(out_jsonl), outdir=str(out_srr_dir))

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
        srrs=test_srrs,                   # TEST SUBSET: two random SRR runs only
        outdir=out_assembly_dir
    )
    # Full pool instead:
    # completed = f.ass.assembly(srrs=pool, outdir=out_assembly_dir)
    
    for d in sorted(out_assembly_dir.iterdir()):
        
        if not d.is_dir():
            continue
        
        keep = {f"{d.name}.fa", f"{d.name}_index"}
        
        for item in d.iterdir():
            
            if item.name in keep:
                continue
            
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
    
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
        srrs=test_srrs,                    # TEST SUBSET: two random SRR runs only
        read_dir=out_srr_dir,
        outdir=out_assembly_dir,
    )
    # Full pool instead:
    # completed = f.ass.map_reads(srrs=pool, read_dir=out_srr_dir,
    #                             outdir=out_assembly_dir)

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
        srrs=test_srrs,                    # TEST SUBSET: two random SRR runs only
        outdir=out_assembly_dir,
        bindir=out_bin_dir,
    )
    # Full pool instead:
    # completed = f.ass.bin_contigs(srrs=pool, outdir=out_assembly_dir,
    #                              bindir=out_bin_dir)

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
    
    gulf_bins = Path("Data/Raw/Bins")
    out_metadata = tmp / "metadata.csv"
    f.metadata.move_bins(out_bin_dir, gulf_bins)
    f.metadata.fetch_metadata(out_jsonl, gulf_bins, out_metadata)