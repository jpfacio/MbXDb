# mxbdb-cli

`mxbdb-cli.py` is a standalone, pandas-only query CLI over the MbXDb entity CSVs
(`Data/Entities/{genes,pah,bins}.csv`). It does not import the pipeline
`functions` package, so it avoids pulling in the full toolchain (habanero,
goatools, biopython) and starts fast.

Run it from the repo root so the relative `Data/Entities/` paths resolve.

```bash
python mxbdb-cli.py <command> [args]
```

## Inputs

| File | Columns used |
| --- | --- |
| `Data/Entities/genes.csv` | `Gene_Tag`, `Contig`, `Start`, `Stop`, `Product`, `Databases`, `Product_Sequence`, `Bin` |
| `Data/Entities/pah.csv` | `Gene_Tag`, `Bin`, `Product`, `Family`, `Subject`, `Identity`, `Evalue`, `Bitscore`, `Qcov` |
| `Data/Entities/bins.csv` | `Bin`, `Sample`, `Project`, `Date`, `Study_ID` |

## Output formats

Every subcommand accepts the same I/O flags:

- `--csv` – print as CSV (no index).
- `--json` – print as JSON records, indented.
- (default) – aligned pretty table via `DataFrame.to_string`.

## Commands

### `genes search <keyword>`

Substring match (case-insensitive, not regex) against `Product` and `Gene_Tag`.

| Option | Description |
| --- | --- |
| `--bin BIN` | Restrict to a single bin. |
| `--limit N` | Show at most N rows. |

Output columns: `Gene_Tag`, `Bin`, `Product`, `Databases`.

```bash
python mxbdb-cli.py genes search dioxygenase --bin GOMC.bin.13102 --limit 10
```

### `genes get <tag>`

Show the full gene record for a `Gene_Tag`, followed by its PAH hits.

```bash
python mxbdb-cli.py genes get GOMC.bin.13102_1547 --json
```

Exits with an error (non-zero) if the tag is not found.

### `pah top`

List PAH hits sorted by `Evalue` (ascending — best hits first). Requires at
least one filter or `-n/--top`, otherwise it errors out.

| Option | Description |
| --- | --- |
| `--family FAMILY` | PAH family (e.g. `nahAc`, `nidA`, `phdF`). |
| `--bin BIN` | Restrict to a single bin. |
| `--min-identity FLOAT` | Minimum `Identity`. |
| `--min-qcov FLOAT` | Minimum `Qcov`. |
| `--max-evalue FLOAT` | Maximum `Evalue`. |
| `--min-bitscore FLOAT` | Minimum `Bitscore`. |
| `-n, --top N` | Show only the N best hits. |

```bash
python mxbdb-cli.py pah top --family nahAc --min-identity 40 -n 5
python mxbdb-cli.py pah top --max-evalue 1e-20 --min-bitscore 200
```

### `pah summary`

Per-`Bin` × `Family` hit counts (`N`), sorted by bin then family.

| Option | Description |
| --- | --- |
| `--family FAMILY` | Restrict to a single family. |
| `--bin BIN` | Restrict to a single bin. |

```bash
python mxbdb-cli.py pah summary --bin GOMC.bin.13102
```

### `bins`

Overview with one row per bin, joining gene and PAH counts onto `bins.csv`.

Output columns: `Bin`, `Sample`, `Study_ID`, `Gene_Count`, `Pah_Hits`,
`Pah_Families`. Counts default to `0` for bins with no genes/hits.

```bash
python mxbdb-cli.py bins --csv
```

## Notes

- The entity CSVs are pipeline outputs and are gitignored; generate them with
  `./mxbdb-build.sh --chen` (or `--all`) before querying.
- Table width/column width are fixed via `pandas` display options
  (`display.width=140`, `display.max_colwidth=60`).
