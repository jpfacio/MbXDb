"""Merge a staging Data tree into the master Data tree.

Each chunked pipeline run writes its own metadata, entities and reports under
a staging Data directory. This module appends those rows to the master CSVs
(deduplicating on the natural key of each table) and copies the raw bins and
Bakta outputs that the chunk produced.
"""

import argparse
import shutil
from pathlib import Path

import pandas as pd


TABLES = [
    ("Entities/bins.csv", ["Bin"]),
    ("Entities/studies.csv", ["Study_ID"]),
    ("Entities/genes.csv", ["Gene_Tag"]),
    ("Entities/annotations.csv", ["Gene_Tag"]),
    ("Entities/pah.csv", ["Gene_Tag", "Family"]),
    ("Entities/pah_annotations.csv", ["Gene_Tag"]),
    ("Reports/pah_summary.csv", ["Bin", "Family"]),
    ("metadata.csv", ["bin"]),
]


def merge_csv(src: Path, dest: Path, keys: list[str]) -> int:
    if not src.exists() or src.stat().st_size == 0:
        return 0

    try:
        new = pd.read_csv(src, encoding="utf-8-sig")
    except pd.errors.EmptyDataError:
        return 0

    if new.empty:
        return 0

    missing = [key for key in keys if key not in new.columns]
    if missing:
        raise ValueError(f"{src}: missing merge key(s) {', '.join(missing)}")

    old = None
    if dest.exists() and dest.stat().st_size > 0:
        try:
            old = pd.read_csv(dest, encoding="utf-8-sig")
        except pd.errors.EmptyDataError:
            old = None

    if old is not None:
        columns = list(old.columns)
        for column in new.columns:
            if column not in columns:
                columns.append(column)
        old = old.reindex(columns=columns)
        new = new.reindex(columns=columns)
        combined = pd.concat([old, new], ignore_index=True)
    else:
        combined = new

    combined = combined.drop_duplicates(subset=keys, keep="last")
    dest.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(dest, index=False)
    return len(combined)


def copy_tree(src: Path, dest: Path) -> int:
    if not src.exists():
        return 0

    copied = 0
    dest.mkdir(parents=True, exist_ok=True)

    for item in sorted(src.iterdir()):
        target = dest / item.name

        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        elif not target.exists():
            shutil.copy2(item, target)

        copied += 1

    return copied


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge a staging Data tree into the master Data tree."
    )
    parser.add_argument("--src", type=Path, required=True,
                        help="staging Data directory produced by a chunk")
    parser.add_argument("--dest", type=Path, required=True,
                        help="master Data directory")
    args = parser.parse_args()

    for relative in ("Raw/Bins", "Raw/Processed"):
        copied = copy_tree(args.src / relative, args.dest / relative)
        if copied:
            print(f"{relative}: {copied:,} items")

    for relative, keys in TABLES:
        rows = merge_csv(args.src / relative, args.dest / relative, keys)
        if rows:
            print(f"{relative}: {rows:,} rows")


if __name__ == "__main__":
    main()
