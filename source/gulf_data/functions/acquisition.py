import subprocess
import json


def get_gcf_info(accession: str) -> list:
    """Fetch the assembled MAG (GCA) summary for a BioProject.

    Used when the BioProject contains assembled metagenomes.

    Args:
        accession (str): BioProject accession (PRJNA...)

    Returns:
        list: One dict per MAG (accession, organism, bioproject, biosample)
    """

    result = subprocess.run(
        ["datasets", "summary", "genome", "accession", accession,
         "--mag", "only", "--as-json-lines"],
        capture_output=True, text=True, check=True
    )

    rows = [json.loads(line) for line in result.stdout.splitlines()]

    subset = [
        {
            "accession": d["accession"],
            "organism": d["organism"]["organism_name"],
            "bioproject": d["assembly_info"]["bioproject_accession"],
            "biosample": d["assembly_info"]["biosample"]["accession"],
        }
        for d in rows
    ]

    return subset


def get_srr_info(accession: str) -> str:
    """Fetch eutils runinfo for an SRR run.

    Used when the BioProject only provides raw reads (no assembled metagenomes).
    Each runinfo row carries SRAStudy, BioProject and BioSample accessions.

    Args:
        accession (str): SRR accession (e.g. "SRR21147283")

    Returns:
        str: The raw runinfo output
    """

    result = subprocess.run(
        ["wget", "-qO-",
         "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
         f"?db=sra&id={accession}&rettype=runinfo&retmode=text"],
        capture_output=True, text=True, check=True
    )

    return result.stdout