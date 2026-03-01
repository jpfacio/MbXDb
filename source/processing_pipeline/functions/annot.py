import subprocess
from pathlib import Path
import pandas as pd
from itertools import islice
import requests
import time

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
    
def create_domain_metadata(genes: str):
    genes = pd.read_csv(genes)
    
    genes_filt = genes[genes['Databases'].notna()]
    domain = genes_filt[['Gene_Tag', 'Databases']].copy()
    
    domain['Databases'] = domain['Databases'].str.split(',')
    
    db_types =[]
    for i in domain['Databases']:
        pref = [s.split(':')[0] for s in i]
        db_types.append(pref)
    
    db_types = list(set(i for sublist in db_types for i in sublist))
    db_types = list(set(s.replace(' ', '') for s in db_types))
    
    for i in db_types:
        domain[i] = None
    
    for idx, db_list in domain['Databases'].items():
        for entry in db_list:
            prefix, value = map(str.strip, entry.split(':', 1))
            if prefix in db_types:
                domain.loc[idx, prefix] = (
                    value if domain.loc[idx, prefix] is None
                    else domain.loc[idx, prefix] + f",{value}"
                )
        
    domain.to_csv("tmp/domain_metadata.csv", index=False)
    

    
    