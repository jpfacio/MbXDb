import pandas as pd
import json

def convert_json_csv(json_file: str) -> None:
    
    """Converts a json input into a csv output
    
    Args:
        json_file (string): Json filepath to be converted
        
        Returns:
            None
    """
    
    with open(json_file, 'r') as f:
        data = json.load(f)

    records = data['data']
    
    df = pd.DataFrame(records)
    df.to_csv("tmp/chen_data.csv")
    
def filtering(csv: str, datadir: str) -> str:
    
    """Filter the input csv and writes a shell script based on its content
    
    Args:
        csv (string): csv filepath
        datadir (string): The target location of files
        
    Returns:
        string: The filepath to the shell script
    """
    df = pd.read_csv(csv)
    total = len(df)
    filtered = df[df["latitude_and_longitude"] != "missing"]
    removed = total - len(filtered)
    
    print(f'Removing {removed} bins without geographical metadata.')
    print('Done!')
    print(f"{len(filtered)} bins left")
    
    filtered["download"] = filtered.apply(
        lambda x: f"{x['download']}{x['assembly_id']}/{x['sample_name']}.fa.gz",
        axis=1
    )
    
    links = filtered["download"].tolist()
    
    wget = [f"wget -P {datadir} -c {link}" for link in links]
    
    with open ("source/chen_data/links.sh", "w") as f:
        f.write('#!/bin/bash\n\n')
        for cmd in wget:
            f.write(cmd + "\n")
            
    return "source/chen_data/links.sh"


    
