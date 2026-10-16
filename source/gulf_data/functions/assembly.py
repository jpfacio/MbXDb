import os
import shutil
import subprocess
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def assembly(srrs: list, 
             outdir: Path = Path("tmp/gulf_tmps/assembly"),
             max_workers: int = 3
) -> int:
    
    srrs = sorted(set(x for x in srrs if x))
    
    to_do = [s for s in srrs if not (outdir / s / f"{s}.fa").exists()]
    
    outdir.mkdir(parents=True, exist_ok=True)
    
    stem = Path("tmp/gulf_tmps/srr")
    paths = []
    for i in to_do:
        
        fwd = i + "_1.fastq"
        rev = i + "_2.fastq"
        
        path_tuple = (stem / fwd, stem / rev)
        
        paths.append(path_tuple)
          
    print(f"Assembling {len(to_do):,} SRRs "
          f"({len(srrs) - len(to_do):,} already present)...\n")
    
    def worker(path):
        
        fwd = path[0]
        rev = path[1]
        srr = fwd.stem.split("_")[0]

        for read in (fwd, rev):
            if not read.exists():
                raise FileNotFoundError(
                    f"reads for {srr} missing ({read.name}) - "
                    "run the download stage first"
                )
        
        target = outdir / srr
        
        if not (target / "done").exists():
            
            if target.exists() and not (target / "options.json").exists():
                shutil.rmtree(target)
            
            subprocess.run(
                ["megahit", "-1", str(fwd), "-2", str(rev), "-o", str(target),
                 "-m", f"{0.9 / max_workers:.2f}",
                 "-t", str(max(1, os.cpu_count() // max_workers)),
                 "--continue"],
                check=True
            )
        
        contigs = target / "final.contigs.fa"
        if not contigs.exists():
            contigs = target / "final_contigs.fa"
        contigs.rename(target / f"{srr}.fa")
        return srr
    
    completed = 0
    if paths:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            
            futures = {
                executor.submit(worker, path): path
                for path in paths
            }
            
            for future in tqdm(
                    as_completed(futures),
                    total=len(paths),
                    desc="megahit",
                    unit="runs"
            ):
                try:
                    future.result()
                    completed += 1
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")
    
    return completed

def build_index(outdir: Path = Path("tmp/gulf_tmps/assembly"),
                max_workers: int = 3
) -> int:
    
    srrs = sorted(d.name for d in outdir.iterdir() if d.is_dir())
    
    contigs = [s for s in srrs if (outdir / s / f"{s}.fa").exists()]
    
    to_do = [s for s in contigs
             if not (outdir / s / f"{s}_index" / f"{s}.rev.2.bt2").exists()]
    
    blocked = [s for s in srrs if s not in set(contigs)]
    
    if blocked:
        print(f"[SKIP] {len(blocked):,} SRRs without contigs "
              f"({', '.join(blocked)})")
    
    print(f"Indexing {len(to_do):,} SRRs "
          f"({len(contigs) - len(to_do):,} already present)...\n")
    
    def worker(srr):
        
        fa = outdir / srr / f"{srr}.fa"
        target = outdir / srr / f"{srr}_index"
        target.mkdir(parents=True, exist_ok=True)
        
        subprocess.run(
            ["bowtie2-build", "--threads",
             str(max(1, os.cpu_count() // max_workers)),
             str(fa), str(target / srr)],
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
                    desc="bowtie2-build",
                    unit="runs"
            ):
                try:
                    future.result()
                    completed += 1
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")
    
    return completed

def map_reads(srrs: list,
              read_dir: Path = Path("tmp/gulf_tmps/srr"),
              outdir: Path = Path("tmp/gulf_tmps/assembly"),
              max_workers: int = 3
) -> int:

    srrs = sorted(set(x for x in srrs if x))

    blocked, to_do = [], []
    for s in srrs:

        fwd = read_dir / f"{s}_1.fastq"
        rev = read_dir / f"{s}_2.fastq"
        bam = outdir / s / f"{s}.bam"
        idx = Path(f"{outdir / s / f'{s}_index' / s}.rev.2.bt2")

        if not idx.exists():
            blocked.append(s)
        elif bam.exists():
            continue
        elif not (fwd.exists() and rev.exists()):
            blocked.append(s)
        else:
            to_do.append(s)

    if blocked:
        print(f"[SKIP] {len(blocked):,} SRRs without index or reads "
              f"({', '.join(blocked)})")

    print(f"Mapping {len(to_do):,} SRRs "
          f"({len(srrs) - len(to_do) - len(blocked):,} already present)...\n")

    threads = max(1, os.cpu_count() // max_workers)

    def worker(srr):

        fwd = read_dir / f"{srr}_1.fastq"
        rev = read_dir / f"{srr}_2.fastq"
        idx = outdir / srr / f"{srr}_index" / srr
        bam = outdir / srr / f"{srr}.bam"
        tmp = outdir / srr / f"{srr}.sorttmp"

        with open(outdir / srr / "bowtie2.log", "w") as err:

            aln = subprocess.Popen(
                ["bowtie2", "-x", str(idx), "-1", str(fwd), "-2", str(rev),
                 "-q", "-p", str(threads)],
                stdout=subprocess.PIPE,
                stderr=err
            )

            sort = subprocess.Popen(
                ["samtools", "sort",
                 "-@", str(threads), "-m", "1G",
                 "-T", str(tmp), "-o", str(bam), "-"],
                stdin=aln.stdout
            )

            aln.stdout.close()
            sort_rc = sort.wait()
            aln_rc = aln.wait()

        if aln_rc != 0 or sort_rc != 0:
            bam.unlink(missing_ok=True)
            raise RuntimeError(
                f"bowtie2 rc={aln_rc} / samtools sort rc={sort_rc}"
            )

        subprocess.run(["samtools", "index", str(bam)], check=True)
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
                    desc="bowtie2",
                    unit="runs"
            ):
                try:
                    future.result()
                    completed += 1
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")

    return completed

def bin_contigs(srrs: list,
                outdir: Path = Path("tmp/gulf_tmps/assembly"),
                bindir: Path = Path("tmp/gulf_tmps/bins"),
                max_workers: int = 3
) -> int:

    srrs = sorted(set(x for x in srrs if x))

    blocked, to_do = [], []
    for s in srrs:

        bam = outdir / s / f"{s}.bam"
        fa = outdir / s / f"{s}.fa"
        target = bindir / s

        if not (bam.exists() and fa.exists()):
            blocked.append(s)
        elif (target / "done").exists():
            continue
        else:
            to_do.append(s)

    if blocked:
        print(f"[SKIP] {len(blocked):,} SRRs without bam or contigs "
              f"({', '.join(blocked)})")

    print(f"Binning {len(to_do):,} SRRs "
          f"({len(srrs) - len(to_do) - len(blocked):,} already present)...\n")

    def worker(srr):

        bam = outdir / srr / f"{srr}.bam"
        fa = outdir / srr / f"{srr}.fa"
        target = bindir / srr

        if target.exists() and any(target.iterdir()):
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)

        depth = target / f"{srr}.depth"

        subprocess.run(
            ["jgi_summarize_bam_contig_depths",
             "--outputDepth", str(depth), str(bam)],
            check=True
        )

        subprocess.run(
            ["metabat2", "-i", str(fa), "-a", str(depth),
             "-o", str(target / srr), "-m", "2500"],
            check=True
        )

        (target / "done").touch()
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
                    desc="metabat2",
                    unit="runs"
            ):
                try:
                    future.result()
                    completed += 1
                except Exception as e:
                    print(f"[ERROR] {futures[future]}: {e}")

    return completed


def _remove_path(path: Path, attempts: int = 3) -> bool:
    """Remove a file or directory, retrying around transient filesystem races.

    MEGAHIT can still be flushing its working directory as it exits, so a
    plain rmtree occasionally hits ENOTEMPTY. Each attempt is separated by a
    short sleep, with a final forced removal as a fallback.

    Returns:
        bool: True when the path no longer exists afterwards
    """
    for _ in range(attempts):
        try:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            return True
        except OSError:
            time.sleep(1)

    try:
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        elif path.exists():
            path.unlink()
    except OSError:
        pass

    return not path.exists()


def clean_assembly_dirs(
    outdir: Path = Path("tmp/gulf_tmps/assembly")
) -> int:
    """Discard MEGAHIT leftovers from each assembly directory.

    Only <srr>.fa and <srr>_index survive; intermediate_contigs/, log and
    options.json are removed. A stubborn item is reported and skipped rather
    than aborting the stage, so the run can continue to mapping/binning.

    Args:
        outdir (Path): root holding one directory per SRR

    Returns:
        int: number of items removed
    """
    outdir = Path(outdir)
    removed = 0

    for d in sorted(outdir.iterdir()):

        if not d.is_dir():
            continue

        keep = {f"{d.name}.fa", f"{d.name}_index"}

        for item in d.iterdir():

            if item.name in keep:
                continue

            if _remove_path(item):
                removed += 1
            else:
                print(f"[WARN] could not remove {item}")

    return removed

