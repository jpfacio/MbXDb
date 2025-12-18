from pathlib import Path
import subprocess
import functions as f

scripts_dir = Path("source/chen_data")
data_dir = Path("Data/Raw")
tmp=Path("tmp")

csv_file = tmp / "chen_data.csv"
json_file =  tmp / "chen_data.json"

subprocess.run([
    "wget", 
    "https://db.cngb.org/maya/api/get_export_data?dataset_id=MDB0000002&table_name=sample",
    "-O",
    str(json_file)
], check=True)
    
f.convert_json_csv(json_file)

f.filtering(csv_file, data_dir / "Bins")

links_script = scripts_dir / "links.sh"

subprocess.run([str(links_script)], check=True)
links_script.unlink()





    
