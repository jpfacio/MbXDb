# MbXDb – Metagenomic biorremediation Xenobiotics Database

v0.0.1

MbXDb is a modular bioinformatics pipeline designed to process metagenome-assembled genomes (MAGs/bins) through quality filtering, functional annotation, and downstream analysis. The pipeline integrates state-of-the-art tools to ensure high-quality, biologically meaningful results.

> **Note**: This project is under active development. Taxonomic and xenobiotic analysis modules are upcoming.

---

## Pipeline Overview

The MXD pipeline executes the following steps sequentially:

1. **Filtering**: Retains only bins with **<1,000 contigs** (as a proxy for complexity/completeness).
2. **Quality Control**: Assesses bin quality using **CheckM2** (completeness, contamination, strain heterogeneity).
3. **Functional Annotation**: Annotates genes and proteins using **Bakta**, providing:
   - COG, Pfam, TIGRFAM annotations
   - rRNA, tRNA, ncRNA detection
   - Signal peptides and transmembrane helices
   - Taxonomic hints via conserved marker genes
4. **(In Progress)** Taxonomic profiling & xenobiotic metabolism potential analysis.

---

## Directory Structure

```
MXD/
├── Data/ # All input and output data
│ ├── Entities/ # Metadata, sample info, entities (reserved)
│ ├── Raw/ # Raw input data
│ │ ├── Bins/ # Input bins
│ │ └── Processed/ # Intermediate filtered bins
│ └── Reports/ # Final analysis reports
│
├── source/ # Analysis code modules
│ ├── ncbi/ # NCBI SRA acquisition/assembly scripts (one run per subproject)
│ ├── processing_pipeline/ # Core pipeline logic 
│ └── test/ # Test subset and validation scripts
│
├── support_files/mining/ncbi/ # NCBI source JSONs (DOI -> SRRs)
├── mxbdb-build.sh # Main entrypoint script 
├── mxbdb-cli.py # Query CLI (Data/Entities CSVs)
└── environment.yml # Conda environment 
```

## Quick Start

### Prerequisites
- Linux
- `conda` (or `mamba`)
- At least 16 GB RAM 

### 1. Clone and set up
```bash
git clone https://github.com/jpfacio/MXD.git
cd MXD
conda env create -f environment.yml
conda activate mxd
```

### 2. Run the pipeline

-Process all NCBI sources
`./mxbdb-build.sh --ncbi`

-Or process a single NCBI source (e.g. gulf)
`./mxbdb-build.sh --ncbi gulf`

-Specify the Bakta database and keep the tmp directory
`./mxbdb-build.sh --ncbi gulf --bakta-db /path/to/db-light --keep-tmps`

-Process the remaining SRRs 10 at a time, removing each batch's tmp data
`./mxbdb-build.sh --ncbi gulf --chunks 10`

-Process a custom study with an inline DOI -> SRR mapping
`./mxbdb-build.sh --ncbi gulf --isolated-run '{"10.xxxx/xxxx": ["SRR123", "SRR456"]}'`

-Or load that mapping from a JSON file
`./mxbdb-build.sh --ncbi gulf --isolated-run /path/to/study.json`

`--isolated-run` always takes a `{doi: [srr]}` JSON object (or a path to one)
and requires exactly one `--ncbi SOURCE`. The SRRs do not need to exist in
`support_files/mining/ncbi/`.

`--chunks` and `--isolated-run` are mutually exclusive and imply `--ncbi` (all
sources unless restricted with `--ncbi SOURCE`). Runs already listed in
`Data/metadata.csv` are skipped, and new results are merged into the existing
`Data` tables. Without either flag the pipeline runs every SRR and then removes
`tmp/` unless `--keep-tmps` is used.

Each batch writes its logs to `log/chunks/<chunk>/` (`run.log`, `bakta_log/`
and `pah_log/`); `log/run.log` indexes the batches. Logs are kept even when the
batch `tmp/` data is removed.
