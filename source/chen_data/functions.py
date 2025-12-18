import pandas as pd
import json

def convert_json_csv(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    records = data['data']
    
    df = pd.DataFrame(records)
    df.to_csv("tmp/chen_data.csv")
    
def filtering(csv, datadir):
    df = pd.read_csv(csv)
    df = df.drop(["sample_type", "environment", "host", "technology", "country", "city",
              "collected_date", "isolation_source", "platform", "altitude", "host_disease"], axis=1)
    df = df[df["latitude_and_longitude"] != "missing"]
    filtered = df.sample(frac=0.05, random_state=5)
    
    filtered.to_csv("tmp/chen_data_5.csv", index=False)
    
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


    
