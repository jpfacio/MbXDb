from pathlib import Path
import subprocess
import functions as f
import pandas as pd

# Defining directories and files

scripts_dir = Path("source/chen_data")
data_dir = Path("Data/Raw/Bins")
tmp=Path("tmp")
csv_file = tmp / "chen_data.csv"
json_file =  tmp / "chen_data.json"

# Download and organizing Data
wget = True

if wget:
    subprocess.run([
        "wget", 
        "https://db.cngb.org/maya/api/get_export_data?dataset_id=MDB0000002&table_name=sample",
        "-O",
        str(json_file)
    ], check=True)
    
    f.convert_json_csv(json_file)

    links_script = f.filtering(csv_file, data_dir)
    
    script_path = Path(links_script)
    script_path.chmod(0o755)

    try:
        subprocess.run([str(script_path)], check=True)
    finally:
        script_path.unlink()
else:
    pass
    
# Organizing the metadata file

chen_data = pd.read_csv("tmp/chen_data.csv")
metadata = pd.read_csv("tmp/metadata.csv")

metadata['bin'] = chen_data['sample_name']
metadata['sample'] = chen_data["assembly_id"]
metadata['project'] = chen_data['project_id']
metadata['id_study'] = "https://doi.org/10.1038/s41586-024-07891-2"
metadata['coord'] = chen_data['latitude_and_longitude']
metadata['date'] = chen_data['collected_date']

metadata.to_csv("tmp/metadata.csv", index=False)









    
