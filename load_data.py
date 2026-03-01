import sqlite3
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

DB_PATH = BASE_DIR / "database" / "mxd.db"
DATA_DIR = BASE_DIR / "Data" / "Entities"


def load_studies(conn):
    path = DATA_DIR / "studies.csv"
    df = pd.read_csv(path)

    df = df.rename(columns={
        "Study_ID": "study_id",
        "Title": "title",
        "Keywords": "keywords",
        "Authors": "authors",
        "Journal": "journal",
        "Year": "year",
        "Affiliations": "affiliations"
    })

    # limpeza
    df["study_id"] = df["study_id"].astype(str).str.strip()

    df = df.drop_duplicates(subset=["study_id"])

    df.to_sql("studies", conn, if_exists="append", index=False)

    print(f"✔ studies carregado ({len(df)} entries)")


def load_bins(conn):
    path = DATA_DIR / "bins.csv"
    df = pd.read_csv(path)

    df = df.rename(columns={
        "Bin": "bin",
        "Sample": "sample",
        "Project": "project",
        "Date": "date",
        "Study_ID": "study_id"
    })

    # limpeza
    df["bin"] = df["bin"].astype(str).str.strip()
    df["study_id"] = df["study_id"].astype(str).str.strip()

    df = df.drop_duplicates(subset=["bin"])

    # 🔍 validação FK (study_id)
    valid_studies = set(pd.read_sql("SELECT study_id FROM studies", conn)["study_id"])
    invalid = df[~df["study_id"].isin(valid_studies)]

    if not invalid.empty:
        print("❌ study_id inválidos em bins:")
        print(invalid["study_id"].unique()[:10])
        raise ValueError("Erro de integridade em bins → studies")

    df.to_sql("bins", conn, if_exists="append", index=False)

    print(f"✔ bins carregado ({len(df)} entries)")


def load_genes(conn):
    path = DATA_DIR / "genes.csv"
    df = pd.read_csv(path)

    df = df.rename(columns={
        "Gene_Tag": "gene_tag",
        "Contig": "contig",
        "Start": "start",
        "Stop": "stop",
        "Product": "product",
        "Databases": "databases",
        "Product Sequence": "product_sequence",
        "Bin": "bin"
    })

    # limpeza
    df["gene_tag"] = df["gene_tag"].astype(str).str.strip()
    df["bin"] = df["bin"].astype(str).str.strip()

    # tipos
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["stop"] = pd.to_numeric(df["stop"], errors="coerce")

    df = df.drop_duplicates(subset=["gene_tag"])

    # 🔍 validação FK (bin)
    valid_bins = set(pd.read_sql("SELECT bin FROM bins", conn)["bin"])
    invalid = df[~df["bin"].isin(valid_bins)]

    if not invalid.empty:
        print("❌ bins inválidos em genes:")
        print(invalid["bin"].unique()[:10])
        raise ValueError("Erro de integridade em genes → bins")

    df.to_sql("genes", conn, if_exists="append", index=False)

    print(f"✔ genes carregado ({len(df)} entries)")


def main():
    conn = sqlite3.connect(DB_PATH)

    conn.execute("PRAGMA foreign_keys = ON;")

    try:
        load_studies(conn)
        load_bins(conn)
        load_genes(conn)

        print("🎉 Banco populado com sucesso")

    except Exception as e:
        print("❌ Error during loading:", e)

    finally:
        conn.close()


if __name__ == "__main__":
    main()