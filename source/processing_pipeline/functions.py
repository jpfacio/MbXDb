import subprocess
from pathlib import Path

def seqkit_summary(bin: str, output):
    bin_path = Path(bin)
    
    subprocess.run([
        "seqkit", "stats",
        "-a", "-T", "-o", str(output),
        str(bin_path)
    ], check=True)
    
def checkm2_analysis(bin: str, output):
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
    
def bakta_analysis(bin: str, db: str, output):
    bin_path = Path(bin)
    
    subprocess.run([
        "bakta", "--db", 
        str(db),
        "--output",
        str(output),
        "--meta", "--threads", "4",
        "--force",
        str(bin)
    ])
    
    
    
    