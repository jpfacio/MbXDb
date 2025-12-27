import pandas as pd

df = pd.read_csv("../../Data/Raw/Processed/GCA_026645595.1_ASM2664559v1_genomic.inference.tsv",
            sep="\t", comment="#", header=None)


df_filtrado = df[
    (df.iloc[:, 1] == 'cds') & 
    (df.iloc[:, 7] != 'hypothetical protein')
]

print(df.head())

