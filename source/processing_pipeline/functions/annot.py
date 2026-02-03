import subprocess
from pathlib import Path

def bakta_analysis(bin: str, db: str, output):
    bin_path = Path(bin)
    
    subprocess.run([
        "bakta", "--db", 
        str(db),
        "--output",
        str(output),
        "--meta", "--threads", "12",
        "--force",
        str(bin)
    ])