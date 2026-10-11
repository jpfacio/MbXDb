#!/bin/bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo"

if [ -f mbxdb.sif ]; then
    echo "mbxdb.sif already exists; nothing to do."
    exit 0
fi

if [ -f mbxdb-docker.tar ]; then
    echo "Building mbxdb.sif from mbxdb-docker.tar (offline)..."
    apptainer build mbxdb.sif docker-archive://mbxdb-docker.tar
else
    echo "Building mbxdb.sif from mbxdb.def (needs internet)..."
    apptainer build mbxdb.sif mbxdb.def
fi

echo "Done: $repo/mbxdb.sif"
