import subprocess
from pathlib import Path
import pandas as pd
from itertools import islice
import requests
import time
import ast

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
    
def map_prot(df_path):
    df = pd.read_csv(df_path)
    
    mapping = (
        df[["UniRef", "Gene_Tag"]]
        .dropna()
        .drop_duplicates(subset="UniRef")
        .set_index("UniRef")["Gene_Tag"]
        .to_dict()
    )

    uniref_list = (
        df["UniRef"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()[:20]
    )

    results_list = []
    url = "https://rest.uniprot.org/uniprotkb/search"

    for i, uniref in enumerate(uniref_list, 1):
        print(f"[{i}/{len(uniref_list)}] Testing UniRef: {uniref}")

        accession = None

        try:
            params = {
                "query": f"uniref_cluster_50:{uniref}",
                "format": "json",
                "size": 5
            }

            r = requests.get(url, params=params)
            r.raise_for_status()

            data = r.json()
            results = data.get("results") or []

            if results:
                accession = results[0].get("primaryAccession")

            if accession is None:
                fallback_id = uniref.replace("UniRef50_", "")
                print(f"  → Fallback to accession: {fallback_id}")

                params = {
                    "query": f"accession:{fallback_id}",
                    "format": "json",
                    "size": 1
                }

                r = requests.get(url, params=params)
                r.raise_for_status()

                data = r.json()
                results = data.get("results") or []

                if results:
                    accession = results[0].get("primaryAccession")

            if accession is None:
                print("  → No results found (both attempts)")
            else:
                print(f"  → Found: {accession}")

        except Exception as e:
            print(f"  → Error: {e}")

        results_list.append({
            "Gene_Tag": mapping.get(uniref),
            "UniRef": uniref,
            "UniProtKB": accession
        })

        time.sleep(0.3) 

    df_out = pd.DataFrame(results_list)
    df_out["UniProtKB"] = df_out["UniProtKB"].astype("object")

    return df_out

def uniprot_info(df_path, batch_size=200):
    df = pd.read_csv(df_path)

    uniprot_list = (
        df["UniProtKB"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if not uniprot_list:
        print("No UniProtKB entries found")
        return pd.DataFrame(columns=[
            "UniProtKB", "Protein_name", "Gene",
            "Organism", "Length", "Reviewed", "GO_terms"
        ])

    url = "https://rest.uniprot.org/uniprotkb/search"
    all_results = []

    for b in range(0, len(uniprot_list), batch_size):
        batch = uniprot_list[b:b + batch_size]

        print(f"\nProcessing batch {b} -> {b + len(batch)}")

        query = " OR ".join([f"accession:{acc}" for acc in batch])

        params = {
            "query": query,
            "fields": ",".join([
                "accession",
                "protein_name",
                "gene_primary",
                "organism_name",
                "length",
                "reviewed",
                "go_id"
            ]),
            "format": "json",
            "size": batch_size
        }

        batch_results = {}

        try:
            r = requests.get(url, params=params)
            r.raise_for_status()
            data = r.json()

            while True:
                entries = data.get("results") or []

                for entry in entries:
                    acc = entry.get("primaryAccession")

                    protein = (
                        entry.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value")
                    )

                    gene = None
                    genes = entry.get("genes") or []
                    if genes:
                        gene = genes[0].get("geneName", {}).get("value")

                    organism = entry.get("organism", {}).get("scientificName")
                    length = entry.get("sequence", {}).get("length")

                    entry_type = entry.get("entryType", "")
                    reviewed = "reviewed" if "reviewed" in entry_type.lower() else "unreviewed"

                    go_terms_list = []

                    for ref in entry.get("uniProtKBCrossReferences", []):
                        if ref.get("database") == "GO":
                            go_terms_list.append(ref.get("id"))

                    go_terms = ";".join(go_terms_list) if go_terms_list else None

                    batch_results[acc] = {
                        "UniProtKB": acc,
                        "Protein_name": protein,
                        "Gene": gene,
                        "Organism": organism,
                        "Length": length,
                        "Reviewed": reviewed,
                        "GO_terms": go_terms
                    }

                next_link = data.get("next")
                if not next_link:
                    break

                r = requests.get(next_link)
                r.raise_for_status()
                data = r.json()

            for acc in batch:
                if acc in batch_results:
                    all_results.append(batch_results[acc])
                else:
                    all_results.append({
                        "UniProtKB": acc,
                        "Protein_name": None,
                        "Gene": None,
                        "Organism": None,
                        "Length": None,
                        "Reviewed": None,
                        "GO_terms": None
                    })

        except Exception as e:
            print(f"Error in batch {b}: {e}")

            for acc in batch:
                all_results.append({
                    "UniProtKB": acc,
                    "Protein_name": None,
                    "Gene": None,
                    "Organism": None,
                    "Length": None,
                    "Reviewed": None,
                    "GO_terms": None
                })

        time.sleep(0.5)

    df_out = pd.DataFrame(all_results)

    expected_cols = [
        "UniProtKB", "Protein_name", "Gene",
        "Organism", "Length", "Reviewed", "GO_terms"
    ]

    df_out = df_out.reindex(columns=expected_cols)

    if not df_out.empty:
        df_out = df_out.astype("object")

    return df_out

def get_go_terms(df_path, batch_size=200):
    df = pd.read_csv(df_path)
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
    
    print("Analysing GO terms")
    
    df = df.dropna(subset=["GO_terms"])
    
    go_dict = {}
    go_ids_list_marked = (
        df.assign(GO_split=df["GO_terms"].str.split(";"))
        .explode("GO_split")
        .assign(combined=lambda x: x["UniProtKB"].str.strip() + "_" + x["GO_split"].str.strip())
        ["combined"]
        .dropna()
        .tolist()
    )
    
    go_ids_list = []
    for i in go_ids_list_marked:
        go_ids_list.append(i.split("_")[1])
    
    
    for i in range(0, len(go_ids_list), batch_size):
        batch = go_ids_list[i:i + batch_size]
        
        print(f"Analyzing GO terms {i} -> {i + len(batch)}")
        
        ids_str = ",".join(batch)
        
        try:
            r = requests.get(f"{base_url}/{ids_str}")
            r.raise_for_status
            data = r.json()
            
            for result in data.get("results", []):
                go_id = result.get("id")
                name = result.get("name")
                
                go_dict[go_id] = name
                
        except Exception as e:
            print(f"Error: {e}")
            
        time.sleep(0.2)

        df_go = pd.DataFrame(
            list(go_dict.items()),
            columns=["go", "Description"]
        )
        
        df_pairs = pd.DataFrame({"pair": go_ids_list_marked})
        df_pairs[["marker", "go"]] = df_pairs["pair"].str.split("_", expand=True)
        
        df_merged = df_pairs.merge(df_go, on="go", how="left")
        
        df_merged = (
            df_merged.groupby("marker", as_index=False)
            .agg({
                "go":list,
                "Description":list
            })
        )
        
        target_df = pd.read_csv("tmp/uniref_to_uniprotkb.csv")
        
        if {"go", "Description"}.issubset(target_df.columns):
            uniprotkb_go = target_df.copy()
        else:
    
    
            uniprotkb_go = target_df.merge(
            df_merged,
            left_on="UniProtKB",
            right_on="marker",
            how="left"
            )
    
    
            if "marker" in uniprotkb_go.columns:
                uniprotkb_go = uniprotkb_go.drop(columns=["marker"])
    
            uniprotkb_go.to_csv("tmp/uniref_to_uniprotkb.csv", index=False)


        genes = pd.read_csv("Data/Entities/genes.csv")
        
        uniprotkb_go = uniprotkb_go[["Gene_Tag", "go", "Description"]]
        
        filt_uniprotkb_go = uniprotkb_go[
            uniprotkb_go["go"].notna() | uniprotkb_go["Description"].notna()
        ].copy()
        
        genes_go = genes.merge(
            filt_uniprotkb_go[["Gene_Tag", "go", "Description"]],
            on="Gene_Tag",
            how="left"
        )
        
        genes_go = genes_go[["Gene_Tag", "Contig", "Start", "Stop",
                             "Product", "Databases", "go", "Description",
                             "Bin", "Product_Sequence"]]
        
        genes_go.to_csv("Data/Entities/genes.csv")
        
def uniparc_info(df_path): 
    df = pd.read_csv(df_path) 
    uniparc_query = df[df["go"].isna()] 
    
    acc_list = uniparc_query["UniProtKB"].dropna().unique().tolist() 
    
    results = [] 
    
    for acc in acc_list: 
        url = f"https://rest.uniprot.org/uniprotkb/{acc}" 
        
        try: 
            r = requests.get(url) 
            r.raise_for_status() 
            data = r.json() 
            
            uniParcId = None 
            if "extraAttributes" in data:
                uniParcId = data["extraAttributes"].get("uniParcId") 
                
                results.append({ 
                                "UniProtKB": acc, 
                                "UniParcId": uniParcId, 
                                "EntryType": data.get("entryType"), 
                                "Raw": data 
                            })
                 
        except Exception as e: 
            print(f"Error in {acc}: {e}") 
            results.append({ 
                            "UniProtKB": acc, 
                            "UniParcId": None, 
                            "EntryType": None, 
                            "Raw": None 
                        }) 
            
        time.sleep(0.2) 
        
    uniparc_info_raw = pd.DataFrame(results) 
    uniparc_info_raw.to_csv("tmp/uniparc_info_raw.csv")
    
    upis = uniparc_info_raw["UniParcId"].dropna().unique().tolist()
    
    results_upi = []
    
    for upi in upis:
        url = f"https://rest.uniprot.org/uniparc/{upi}"
        
        try:
            r = requests.get(url)
            r.raise_for_status()
            data = r.json()
            
            results_upi.append({
                "UniParcId": upi,
                "Raw": data
            })
            
        except Exception as e:
            print(f"Error in {upi}: {e}")
            results_upi.append({
                "UniParcId": upi,
                "Raw": None
            })
    
        time.sleep(0.2)
        
        uniparc_info = pd.DataFrame(results_upi)
        
        interpro_list = []
        
        for raw in uniparc_info["Raw"]:
            if pd.isna(raw):
                interpro_list.append(None)
                continue
            
            if isinstance(raw, str):
                raw = ast.literal_eval(raw)
                
            ids = []
            
            for i in raw.get("sequenceFeatures", []):
                interpro = i.get("interproGroup", {})
                ip_id = interpro.get("id")
                
                if ip_id:
                    ids.append(ip_id)
                    
            ids = ";".join(sorted(set(ids))) if ids else None
            
            interpro_list.append(ids)
            
        uniparc_info["InterPro"] = interpro_list
        uniparc_info = uniparc_info[["UniParcId", "InterPro", "Raw"]]
        uniparc_info.to_csv("tmp/uniparc_info.csv")
        
def interpro_to_go(df_path):
    df = pd.read_csv("tmp/uniparc_info.csv")
    
    
    
        
        
        
           
            
                
                