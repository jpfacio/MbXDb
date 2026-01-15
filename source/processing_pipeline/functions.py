import subprocess
from pathlib import Path
import pandas as pd

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
        "--meta", "--threads", "12",
        "--force",
        str(bin)
    ])
    
def create_genes_ent():
    df = pd.read_csv("Data/Raw/Processed/GCA_026645595.1_ASM2664559v1_genomic.tsv",
            sep="\t")

    gene_list = []
    n_row = len(df['SequenceId'])
    str_base = "GCA_026645595.1"
    tags = [f"{str_base}_{i:04d}" for i in range(1, n_row + 1)]

    genes["Gene_Tag"] = tags
    genes["Contig"] = df["SequenceId"]
    genes["Start"] = df["Start"]
    genes["End"] = df["Stop"]
    genes["Function"] = df["Product"]
    genes["Bin"] = str_base
    genes["Sample_ID"] = "SAMN31654591"

    genes.to_csv("Data/Entities/genes.csv")
    



    
    
    
