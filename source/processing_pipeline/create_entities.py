import pandas as pd

df = pd.read_csv("../../Data/Raw/Processed/GCA_026645595.1_ASM2664559v1_genomic.tsv",
            sep="\t")

gene_list = []
n_row = len(df['SequenceId'])
str_base = "GCA_026645595.1"
sample_id = "SAMN31654591"
project_id = "PRJNA899357"
doi = "https://doi.org/10.1016/j.envpol.2022.120772"
tags = [f"{str_base}_{i:04d}" for i in range(1, n_row + 1)]

genes = pd.DataFrame(columns=["Gene_Tag", "Contig", "Start", "End", 
                              "Function", "Bin", "Sample_ID"])

genes["Gene_Tag"] = tags
genes["Contig"] = df["SequenceId"]
genes["Start"] = df["Start"]
genes["End"] = df["Stop"]
genes["Function"] = df["Product"]
genes["Bin"] = str_base
genes["Sample_ID"] = "SAMN31654591"

genes.to_csv("../../Data/Entities/genes.csv")

bins = pd.DataFrame([{"Bins": str_base,
                              "Sample_ID": sample_id,
                              "Project_ID": project_id,
                              "Study_ID": doi.strip()
                              }])

bins.to_csv("../../Data/Entities/bins.csv", index=False)

title = "Microbiome enrichment from contaminated marine sediments unveils novel bacterial strains for petroleum hydrocarbon and heavy metal bioremediation"
keywords = "Polycyclic aromatic hydrocarbonsHeavy metalsBioremediationNext-generation sequencingMarine biotechnology"
authors = "F Dell'Anno, LJ van Zyl, M Trindade, E Buschi, A Cannavacciuolo,"
year = "2023"
country = "Italy"

studies = pd.DataFrame([{"Study_ID": doi,
                         "Title": title, 
                         "Keywords": keywords,
                         "Authors": authors,
                         "Year": year,
                         "Country": country}])

studies.to_csv("../../Data/Entities/studies.csv", index=False)

coord = "40.73 N 14.47 E"
region = "Tyrrhenian Sea"

env_characteristics = pd.DataFrame([{"Coordinates": coord,
                                     "Region": region,
                                     "Temperature": "",
                                     "Salinity": "",
                                     "Depth (m)": 0.2,
                                     "Sample_ID": sample_id,
                                     "Study_ID": doi}])

env_characteristics.to_csv("../../Data/Entities/env_characteristics.csv")









