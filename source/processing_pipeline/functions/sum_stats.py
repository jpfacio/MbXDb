import subprocess
from pathlib import Path

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
        f"seqkit stats -a -T -o {out_file} {data_dir}/*.fa.gz",
        shell=True,
        check=True
    )
    
def checkm2_analysis(bin: str, output: str) -> None:
    
    """Run CheckM2  Analysis
    
    Args:
        bin (string): Files containing the sequence
        output (string): Output filepath
        
    Returns:
        None
        
    """
    
    bin_path = Path(bin)
    
    subprocess.run([
        "checkm2", "predict", "--input",
        str(bin_path),
        "--output-directory",
        str(output),
        "--threads 4",
        "-x fa",
        "--force"        
    ], check=True)