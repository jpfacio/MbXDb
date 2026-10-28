#!/bin/bash

run_ncbi=false
ncbi_source="all"
keep_tmps=false
bakta_db=""
chunks=""
chunks_set=false
isolated_run=""
isolated_set=false

ncbi_dir="support_files/mining/ncbi"

usage() {
    cat << 'EOF'
Usage: ./mxbdb-build.sh [--ncbi [SOURCE[,SOURCE...]|all]]
                        [--chunks N | --isolated-run VALUE]
                        [--bakta-db DIR] [--keep-tmps]

  --ncbi [SOURCE]  Process NCBI source(s). SOURCE matches a JSON file in
                   support_files/mining/ncbi/ and a run_SOURCE.py script in
                   source/ncbi/ (e.g. gulf). Defaults to all sources.
  --chunks N       Process the remaining SRRs N at a time, removing the tmp
                   directory after each batch.
  --isolated-run V Run a single batch from a {doi: [srr]} JSON object or a
                   path to a JSON file. Requires one --ncbi SOURCE.
  --bakta-db DIR   Path to the Bakta db-light directory. Defaults to
                   ../db-light.
  --keep-tmps      Keep the tmp directory after the run (removed by default).
EOF
}

# Defining arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            usage
            exit 0
            ;;
        --ncbi)
            run_ncbi=true
            shift
            if [[ $# -gt 0 && "$1" != -* ]]; then
                ncbi_source="$1"
                shift
            fi
            ;;
        --ncbi=*)
            run_ncbi=true
            ncbi_source="${1#*=}"
            shift
            ;;
        --bakta-db)
            shift
            if [[ $# -eq 0 || "$1" == -* ]]; then
                echo "Error: --bakta-db needs a directory." >&2
                exit 1
            fi
            bakta_db="$1"
            shift
            ;;
        --bakta-db=*)
            bakta_db="${1#*=}"
            shift
            ;;
        --chunks)
            chunks_set=true
            shift
            if [[ $# -eq 0 || "$1" == -* ]]; then
                echo "Error: --chunks needs a number." >&2
                exit 1
            fi
            chunks="$1"
            shift
            ;;
        --chunks=*)
            chunks_set=true
            chunks="${1#*=}"
            shift
            ;;
        --isolated-run)
            isolated_set=true
            shift
            if [[ $# -eq 0 || "$1" == -* ]]; then
                echo "Error: --isolated-run needs a {doi: [srr]} JSON object or JSON file." >&2
                exit 1
            fi
            isolated_run="$1"
            shift
            ;;
        --isolated-run=*)
            isolated_set=true
            isolated_run="${1#*=}"
            shift
            ;;
        --keep-tmps)
            keep_tmps=true
            shift
            ;;
        *)
            echo "Invalid argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if $chunks_set && $isolated_set; then
    echo "Error: --chunks and --isolated-run cannot be combined." >&2
    exit 1
fi

if $chunks_set; then
    if [[ ! "$chunks" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: --chunks expects a positive integer." >&2
        exit 1
    fi
    run_ncbi=true
fi

if $isolated_set; then
    if [[ -z "${isolated_run// /}" ]]; then
        echo "Error: --isolated-run needs a {doi: [srr]} JSON object or JSON file." >&2
        exit 1
    fi
    run_ncbi=true
fi

if [[ -n "$bakta_db" && ! -d "$bakta_db" ]]; then
    echo "Error: Bakta db directory not found: $bakta_db" >&2
    exit 1
fi

ncbi_sources=()
if $run_ncbi; then
    if [[ "$ncbi_source" == "all" ]]; then
        shopt -s nullglob
        for json_file in "$ncbi_dir"/*.json; do
            ncbi_sources+=("$(basename "$json_file" .json)")
        done
        shopt -u nullglob

        if [[ ${#ncbi_sources[@]} -eq 0 ]]; then
            echo "Error: no NCBI source JSONs found in $ncbi_dir" >&2
            exit 1
        fi
    else
        IFS=',' read -r -a ncbi_sources <<< "$ncbi_source"
    fi

    if [[ ${#ncbi_sources[@]} -eq 0 ]]; then
        echo "Error: no NCBI sources selected" >&2
        exit 1
    fi

    for source in "${ncbi_sources[@]}"; do
        if [[ ! -f "$ncbi_dir/$source.json" ]]; then
            echo "Error: missing source JSON: $ncbi_dir/$source.json" >&2
            exit 1
        fi
        if [[ ! -f "source/ncbi/run_$source.py" ]]; then
            echo "Error: missing run script: source/ncbi/run_$source.py" >&2
            exit 1
        fi
    done

    if $isolated_set && [[ ${#ncbi_sources[@]} -ne 1 ]]; then
        echo "Error: --isolated-run requires exactly one --ncbi SOURCE." >&2
        exit 1
    fi
fi

echo "Starting"
echo "Creating Data Directories"

data_dir="Data"
tmp_dir="tmp"
log_dir="log"

sub_data_dirs=("$data_dir/Entities" "$data_dir/Raw" "$data_dir/Raw/Bins" "$data_dir/Raw/Processed" "$data_dir/Reports") 

if [ -d "$data_dir" ]; then
    echo "Directory '$data_dir' already exists."
    read -p "Replace? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Removing existing directory..."
            if rm -rf "$data_dir"; then
                echo "Data directories removed successfully."
            
            elif rm -rf "$tmp_dir"; then
                echo "Tmp directory removed sucessfully."

            elif rm -rf "$log_dir"; then
                echo "Log directory removed sucessfully."

            else
                echo "Error removing the directories."
                exit 1
            fi
            ;;
        *)
            echo "Keeping existing directory. Skipping creation."
            skip_data=true
            ;;
    esac
fi

if [ -z "$skip_data" ]; then
    for dir in "${sub_data_dirs[@]}"; do
        if mkdir -p "$dir"; then
            echo "Created: $dir"
        else
            echo "Failed to create: $dir"
            exit 1
        fi
    done

    if mkdir -p "$tmp_dir"; then
        echo "Tmp directory created"

        cat > "$data_dir"/metadata.csv <<EOF
bin,srr,sample,project,study_id,date,coord,depth
EOF
    else
        echo "Failed to create: $tmp_dir"
    fi
fi

if mkdir -p "$log_dir"; then
    echo "Log directory created"

    cat > "$log_dir"/run.log <<EOF
##########     LOG FILE     ##########

This file records the duration of each checkpoint, memory and space usage, it is intended for developers,
feel free to delete it when the database is built on your machine.

Run: $(date)

Process: $0


EOF

else
    echo "Failed to create: $log_dir"
fi

if $chunks_set || $isolated_set; then
    if $isolated_set; then
        mode="isolated"
        value="$isolated_run"
    else
        mode="chunks"
        value="$chunks"
    fi

    if $chunks_set && ! $keep_tmps; then
        echo "Removing stale tmp directory"
        rm -rf "$tmp_dir"
    fi
    mkdir -p "$tmp_dir"

    sources_csv=$(IFS=,; echo "${ncbi_sources[*]}")

    orch_args=(
        --mode "$mode"
        --value "$value"
        --sources "$sources_csv"
        --json-dir "$ncbi_dir"
        --data-dir "$data_dir"
        --tmp-root "$tmp_dir"
        --log-dir "$log_dir"
    )
    if $keep_tmps; then
        orch_args+=(--keep-tmps)
    fi
    if [[ -n "$bakta_db" ]]; then
        orch_args+=(--bakta-db "$bakta_db")
    fi

    if ! python "source/ncbi/orchestrate.py" "${orch_args[@]}"; then
        echo "Error: chunked NCBI stage failed" >&2
        exit 1
    fi
else
    if $run_ncbi; then
        for source in "${ncbi_sources[@]}"; do
            if ! python3 "source/ncbi/run_$source.py" "$ncbi_dir/$source.json"; then
                echo "Error: ncbi/$source data stage failed" >&2
                exit 1
            fi
        done
    fi

    pp_args=()
    if [[ -n "$bakta_db" ]]; then
        pp_args+=(--bakta-db "$bakta_db")
    fi

    if ! python "source/processing_pipeline/run_pp.py" "${pp_args[@]}"; then
        echo "Error: processing pipeline failed" >&2
        exit 1
    fi
fi

if $keep_tmps; then
    echo "Keeping tmp directory"
else
    echo "Removing tmp directory"
    rm -rf "$tmp_dir"
fi










