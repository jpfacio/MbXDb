#!/bin/bash
#SBATCH --job-name=mbxdb
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
# TODO: set your cluster's partition/account, e.g.
#SBATCH --partition=CPU
#SBATCH --account=9789751

# Usage:
#   sbatch run-slurm.sh --chen
#   sbatch run-slurm.sh --gulf
#   sbatch run-slurm.sh --all
# (the flag after the script name is forwarded to mxbdb-build.sh)

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

ARGS="${*:---chen}"

echo "Job $SLURM_JOB_ID on $(hostname), args: $ARGS"
echo "Start: $(date)"

./run-container.sh $ARGS

echo "Done: $(date)"
