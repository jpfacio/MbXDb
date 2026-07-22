import subprocess
from pathlib import Path
import pandas as pd

def seqkit_summary(data_dir: str, output: str) -> None:
    
    """Run seqkit analysis on all .fa.gz files in a directory
    
    Args:
        data_dir (string): The data filepath containing the sequences
        output (string): Summary files
        
    Returns:
        None
    """
    
    out_file = Path(output) / "summary_stats.tsv"
    
    subprocess.run(
        f"seqkit stats -a -T -o {out_file} -j 12 {data_dir}/*.fa.gz",
        shell=True,
        check=True
    )
    

    
    
    
