from pathlib import Path
import functions as f
import pandas as pd
import requests
import time
import json

seqkit_key = False
bakta_key=False
checkm_key=False
uniprotkb_key = False
basic_info_key = False

data_dir=Path("Data/Raw/Bins")
tmp=Path("tmp")

print("Generating summary statistics")

if seqkit_key: 
    for bin in data_dir.glob("*"):
        output = tmp / f"{bin.stem}_stats.tsv"
        print(f"Processing: {bin.name}")
        f.st.seqkit_summary(str(bin), str(output))

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

    print("Starting CheckM2 analysis")
else:
    pass

if checkm_key:
    for bin in data_dir.glob("*"):
        output = tmp
        print(f"Processing {bin.name}")
        f.st.checkm2_analysis(str(bin), str(output))
else:
    pass
    
print("Starting protein annotation (Bakta)")

if bakta_key:
    for bin in data_dir.glob("*"):
        output = "Data/Raw/Processed"
        print(f"Annotating: {bin.name}")
        f.annot.bakta_analysis(str(bin), str("/home/joao_facio/joao_facio_nfs/db/db-light"), str(output))
else:
    pass

print("Processing protein annotation")

data_annot = Path("Data/Raw/Processed")

tsv_annot = []
for i in data_annot.glob("*.tsv"):
    name = i.name
    if "hypotheticals" not in name and "inference" not in name:
        tsv_annot.append(i)

genes_bins = []
for i in tsv_annot:
    genes_temp = f.ent.create_genes_ent(i)
    if not genes_temp.empty:
        genes_bins.append(genes_temp)
        
if genes_bins:
    genes_ent = pd.concat(genes_bins, ignore_index=True)
    genes_ent.to_csv("Data/Entities/genes.csv", index=False)
    
        
metadata = "tmp/metadata.csv"

f.ent.create_bins_ent(metadata)

bins = pd.read_csv("Data/Entities/bins.csv")

doi_list = bins["Study_ID"].unique().tolist()

doi_dfs = []
for doi in doi_list:
    df = f.ent.create_studies_ent(doi)
    doi_dfs.append(df)

final_doi_df = pd.concat(doi_dfs, ignore_index=True)

final_doi_df.to_csv("Data/Entities/studies.csv", index=False)

f.annot.create_domain_metadata("Data/Entities/genes.csv")

if uniprotkb_key:
    uniref_to_uniprotkb = f.annot.map_prot("tmp/domain_metadata.csv")
    uniref_to_uniprotkb.to_csv("tmp/uniref_to_uniprotkb.csv")
else:
    pass

if basic_info_key:
    uniprotkb_info = f.annot.uniprot_info("tmp/uniref_to_uniprotkb.csv")
    uniprotkb_info.to_csv("tmp/uniprotkb_info.csv")
else:
    pass

go_terms_info = f.annot.get_go_terms("tmp/uniprotkb_info.csv")

f.annot.uniparc_info("tmp/uniref_to_uniprotkb.csv")




   



    

    




