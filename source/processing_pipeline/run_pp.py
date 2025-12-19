from pathlib import Path
import subprocess
import functions as f
import pandas as pd

data_dir=Path("Data/Raw/Bins")
scripts_dir=Path("source/test")
tmp=Path("tmp")

print("Generating summary statistics")

for bin in data_dir.glob("*"):
    output = tmp / f"{bin.stem}_stats.tsv"
    print(f"Processing: {bin.name}")
    f.seqkit_summary(str(bin), str(output))

seqkit_tsv = next(tmp.glob("*.tsv"))

seqkit_df = pd.read_csv(seqkit_tsv, sep="\t")
seqkit_df = seqkit_df[seqkit_df["num_seqs"] > 1000]

remove_count = 0
for i in seqkit_df['file']:
    file_path = data_dir / i 
    if file_path.exists():
        file_path.unlink()
        remove_count += 1
    else:
        print(f"File not found: {file_path}")
        
print(f"{remove_count} files removed")

#print("Starting CheckM2 analysis")


#####RESOLVER DEPOIS#######


#for bin in data_dir.glob("*"):
    #output = tmp
    #print(f"Processing {bin.name}")
    #f.checkm2_analysis(str(bin), str(output))
    
print("Starting protein annotation (Bakta)")

for bin in data_dir.glob("*"):
    output = "Data/Raw/Processed"
    print(f"Annotating: {bin.name}")
    f.bakta_analysis(str(bin), str("/temporario2/15402906/LBUEL-H11-Exploration/db"), str(output))


   



    

    




