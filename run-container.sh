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
if [ -d "$repo/../db-light" ]; then
    binds="$binds,$repo/../db-light"
else
    echo "Warning: $repo/../db-light not found (Bakta DB); binding repo only." >&2
fi

cd "$repo"

exec apptainer exec \
    --bind "$binds" \
    "$sif" \
    bash -c 'export PATH=/opt/conda/envs/mbxdb-env/bin:$PATH; exec ./mxbdb-build.sh "$@"' _ "$@"
