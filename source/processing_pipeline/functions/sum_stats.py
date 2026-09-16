import subprocess
from pathlib import Path
from io import StringIO
import pandas as pd
import os
import time

def seqkit_summary(data_dir: str, output: str) -> pd.DataFrame:
    
    """Run seqkit analysis on all .fa.gz files in a directory
    
    Args:
        data_dir (string): The data filepath containing the sequences
        output (string): Summary files
        
    Returns:
        Summary statistics in dataframe
    """
    
    out_file = Path(output) / "summary_stats.tsv"
    
    subprocess.run(
        f"seqkit stats -a -T -o {out_file} -j 12 {data_dir}/*.fa.gz",
        shell=True,
        check=True
    )
    
    df = pd.read_csv(out_file, sep='\t')
    
    return df
    
def seq_filter(df: pd.DataFrame) -> None:
    
    """Filter a dataframe removing rows where 'num_seqs" values are bigger than 1000
    
    Args:
        df (pd.DataFrame): The data frame to be processed
        
    Returns:
        None
    """
    
    to_remove = df.loc[df['num_seqs'] > 1000, 'file'].to_list()
    
    removed = 0
    
    for file in to_remove:
        path = Path(file)
        if path.exists():
            path.unlink()
            removed += 1
            
    print(f"{removed} files removed")
    
def seq_remove_500(files_dir: Path):
    
    """Check fa.gz files and apply seqkit to remove contigs below 500bp length, uncompressing the file
    in the end.
    
    Args:
        files_dir (Path): Directory containing .fa.gz files
    Returns:
        None
    """
    
    for fa in files_dir.glob("*.fa.gz"):
        
        original_contigs = int(
            subprocess.check_output(
                ['seqkit', 'stats', '-T', str(fa)],
                text=True
            ).strip().split("\n")[1].split("\t")[3]
        )
        
        tmp = fa.with_suffix(".tmp.fa.gz")
        
        with open(tmp, 'wb') as f:
            subprocess.run(
                ['seqkit', 'seq', '-m', '500', str(fa)],
                stdout=f,
                stderr=subprocess.DEVNULL,
                check=True
            )
            
        processed_contigs = int(
            subprocess.check_output(
                ['seqkit', 'stats', '-T', str(tmp)],
                text=True
            ).strip().split("\n")[1].split("\t")[3]
        )
        
        removed = original_contigs - processed_contigs
        
        Path.replace(tmp, fa)
        
        fa.rename(fa.with_suffix(''))
        
        print(f'{fa.name}: {removed} contigs from {original_contigs}.')

def run_checkm(
    files_dir: Path,
    output_dir: Path,
    log: Path,
    threads: int = 8,
    batch_size: int = 1000,
    out_root: Path = None
) -> pd.DataFrame:
    
    """Runs CheckM in memory-bounded, resumable batches and merges the QA results
    
    Bins are symlinked into subdirectories of at most `batch_size` genomes each;
    lineage_wf runs once per batch and per-batch QA tables are concatenated.
    A batch whose QA table already exists is skipped, so interrupted runs resume
    without recomputing finished batches.
    
    Args:
        files_dir (Path): Path containing the .fa files to be processed
        output_dir (Path): Directory where the CheckM output will be sent
        log (Path): Path of the log directory
        threads (int): Number of threads to be passed to each CheckM batch
        batch_size (int): Maximum number of bins processed per batch
        out_root (Path): Optional override for the batching work directory

    Returns:
        pd.DataFrame: Dataframe containing the merged CheckM results
    """
    
    if out_root is None:
        out_root = output_dir / 'checkm_batches'
    
    out_root.mkdir(parents=True, exist_ok=True)
    
    log_file = log / 'checkm.log'
    
    bins = sorted(files_dir.glob('*.fa'))
    
    batches = [bins[i:i + batch_size] for i in range(0, len(bins), batch_size)]
    
    frames = []
    
    with open(log_file, 'a') as log_handle:
        
        for index, batch in enumerate(batches, start=1):
            
            name = f'batch_{index:03d}'
            batch_dir = out_root / f'{name}_in'
            out_batch = out_root / name
            qa_file = out_root / f'{name}_qa.tsv'
            
            if qa_file.exists():
                print(f'{name}: completed QA table found, skipping.')
                
                log_handle.write(f'=== {name}: completed QA table found, skipping ===\n')
                log_handle.flush()
            
            else:
                log_handle.write(f'=== {name} START ({len(batch)} bins) {time.ctime()} ===\n')
                log_handle.flush()
                
                batch_dir.mkdir(parents=True, exist_ok=True)
                
                for fa in batch:
                    link = batch_dir / fa.name
                    if not link.exists():
                        os.symlink(os.path.relpath(fa.resolve(), batch_dir.resolve()), link)
                
                subprocess.run(
                    [
                        'checkm', 'lineage_wf', '-t',
                        str(threads),
                        '-x', 'fa',
                        str(batch_dir),
                        str(out_batch)
                    ],
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                    text=True
                )
                
                with open(qa_file, 'w') as qa_out:
                    subprocess.run(
                        [
                            "checkm", "qa", "--tab_table", "-o", "2",
                            str(out_batch / "lineage.ms"),
                            str(out_batch)
                        ],
                        stdout=qa_out,
                        stderr=log_handle,
                        check=True,
                        text=True
                    )
                
                log_handle.write(f'=== {name} END {time.ctime()} ===\n')
                log_handle.flush()
            
            with open(qa_file) as fh:
                table_text = ''.join(line for line in fh if not line.startswith('['))
            
            frames.append(pd.read_csv(StringIO(table_text), sep='\t'))
    
    final_table = "tmp/checkm_results.tsv"
    
    df = pd.concat(frames, ignore_index=True)
    
    df.to_csv(final_table, sep='\t', index=False)
    
    print("Done!")
    
    return df

def final_processing(seqkit: pd.DataFrame, checkm: pd.DataFrame, tmp_file: Path) -> None:
    
    seqkit = seqkit[
        [
            "file",
            "num_seqs",
            "N50", 
            "GC(%)"
        ]
    ]
    
    seqkit.loc[:, 'file'] = seqkit['file'].str.removesuffix('.gz')
    seqkit = seqkit.loc[seqkit['num_seqs'] <= 1000].copy()
    
    seqkit = seqkit.assign(
        _key=seqkit['file'].apply(lambda p: Path(p).name).str.removesuffix('.fa')
    ).sort_values('_key').reset_index(drop=True)
    
    checkm = checkm.sort_values(checkm.columns[0]).reset_index(drop=True)
    
    bin_ids = checkm[checkm.columns[0]].astype(str)
    
    checkm = checkm[
        [
            "Completeness",
            "Contamination",
            "Strain heterogeneity"
        ]
    ]
    
    if len(bin_ids) != len(seqkit):
        raise ValueError(
            f'CheckM/seqkit row mismatch: {len(bin_ids)} CheckM bins vs '
            f'{len(seqkit)} seqkit rows.'
        )
    
    if not (bin_ids.to_numpy() == seqkit['_key'].to_numpy()).all():
        raise ValueError('CheckM and seqkit bin order mismatch after sorting.')
    
    df = pd.concat([seqkit.drop(columns='_key').reset_index(drop=True), 
                    checkm.reset_index(drop=True)], 
                   axis=1)
    
    df['MIMAG'] = 'Low'
    
    df.loc[
        (df["Completeness"] >= 50) &
        (df["Contamination"] <= 10), "MIMAG"
    ] = "Medium"
    
    df.loc[
        (df["Completeness"] >= 90) &
        (df["Contamination"] <= 5), "MIMAG"
    ] = "High"
    
    low = df.loc[df['MIMAG'] == 'Low', 'file'].to_list()
    
    for file in low:
        subprocess.run(['rm', '-f', file],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        
    df = df.loc[df['MIMAG'] != "Low"]
        
    print(f"{len(low)} low quality MAGs removed.")
    
    table_path = tmp_file / "qc_metrics.tsv"
    
    df.to_csv(table_path, sep="\t", index=False)
    
    subprocess.run(["rm", "tmp/summary_stats.tsv", "tmp/checkm_results.tsv"], check=True)
    


    

    
    
    
