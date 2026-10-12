#!/bin/bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo"

# Locate the container runtime (apptainer preferred, singularity fallback).
runtime=""
for r in "${CONTAINER_RUNTIME:-}" apptainer singularity \
         /usr/bin/apptainer /usr/bin/singularity \
         /usr/local/bin/apptainer /usr/local/bin/singularity; do
    if [ -n "$r" ] && command -v "$r" >/dev/null 2>&1; then
        runtime="$r"
        break
    fi
done

if [ -z "$runtime" ]; then
    echo "Error: neither apptainer nor singularity found." >&2
    echo "Try: 'module load apptainer', or set CONTAINER_RUNTIME=/path/to/runtime." >&2
    exit 1
fi
echo "Runtime: $runtime"

if [ -f mbxdb.sif ]; then
    echo "mbxdb.sif already exists; nothing to do."
    exit 0
fi

if [ -f mbxdb-docker.tar ]; then
    echo "Building mbxdb.sif from mbxdb-docker.tar (offline)..."
    "$runtime" build --fakeroot mbxdb.sif docker-archive://mbxdb-docker.tar
else
    echo "Building mbxdb.sif from mbxdb.def (needs internet)..."
    "$runtime" build --fakeroot mbxdb.sif mbxdb.def
fi

echo "Done: $repo/mbxdb.sif"
