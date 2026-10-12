#!/bin/bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sif="$repo/mbxdb.sif"

if [ ! -f "$sif" ]; then
    echo "Image not found: $sif" >&2
    echo "Build it first:  ./build-container.sh" >&2
    exit 1
fi

binds="$repo"

# Bakta database (gitignored, repo-sibling)
if [ -d "$repo/../db-light" ]; then
    binds="$binds,$repo/../db-light"
else
    echo "Warning: $repo/../db-light not found (Bakta DB); Bakta will fail." >&2
fi

# CheckM1 reference data (~1.4 GB). Not shipped in the image; bind a local copy
# and point CheckM at it. Override with CHECKM_DATA=/path/to/checkm_data.
checkm_src="${CHECKM_DATA:-}"
if [ -z "$checkm_src" ]; then
    for cand in \
        "$HOME/micromamba/envs/mbxdb-env/checkm_data" \
        "${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/mbxdb-env/checkm_data" \
        "$HOME/miniconda3/envs/mbxdb-env/checkm_data" \
        "$repo/../checkm_data"; do
        if [ -f "$cand/taxon_marker_sets.tsv" ]; then
            checkm_src="$cand"
            break
        fi
    done
fi

exec_args=(exec --bind "$binds")
if [ -n "$checkm_src" ]; then
    binds="$binds,$checkm_src:/checkm_data"
    exec_args=(exec --bind "$binds" --env CHECKM_DATA_PATH=/checkm_data)
    echo "CheckM data: $checkm_src -> /checkm_data"
else
    echo "Warning: CheckM reference data not found; set CHECKM_DATA=/path/to/checkm_data (QC stage will fail otherwise)." >&2
fi

cd "$repo"

exec apptainer "${exec_args[@]}" \
    "$sif" \
    bash -c 'export PATH=/opt/conda/envs/mbxdb-env/bin:$PATH; exec ./mxbdb-build.sh "$@"' _ "$@"
