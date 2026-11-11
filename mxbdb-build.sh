#!/bin/bash

run_ncbi=false
run_ncbi_stage=true
run_pp_stage=true
only_ncbi=false
only_pp=false
keep_tmps=false
chunks=""
ncbi_inputs=()

# Defining arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)
            cat << 'EOF'
Usage: ./mxbdb-build.sh --ncbi [name ...] [--only-ncbi | --only-pp]
                      [--chunks N] [--keep-tmps]

  --ncbi [name ...]   Build from the NCBI SRR inputs in support_files/ncbi.
                      Pass one or more json names (e.g. --ncbi gulf arctic).
                      With no names, every json in support_files/ncbi is used.
                      Runs both the NCBI stage and the processing pipeline.
  --only-ncbi         With --ncbi, run only the NCBI stage.
  --only-pp           With --ncbi, run only the processing pipeline.
  --chunks N          Process N SRRs at a time through the NCBI stage, clearing
                      tmp between chunks. Ignored when the NCBI stage is skipped.
  --keep-tmps         Keep the tmp directory after the pipeline finishes.
                      Without it, tmp is removed at the end of the run.
EOF
            exit 0
            ;;
        --only-ncbi)
            only_ncbi=true
            shift
            ;;
        --only-pp)
            only_pp=true
            shift
            ;;
        --chunks)
            if [[ ! "$2" =~ ^[0-9]+$ ]] || [ "$2" -lt 1 ]; then
                echo "Error: --chunks needs a positive integer (got '$2')" >&2
                exit 1
            fi
            chunks="$2"
            shift 2
            ;;
        --keep-tmps)
            keep_tmps=true
            shift
            ;;
        --ncbi)
            run_ncbi=true
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                ncbi_inputs+=("$1")
                shift
            done
            ;;
        *)
            echo "Invalid argument: $1" >&2
            exit 1
            ;;
    esac
done

if ! $run_ncbi; then
    echo "Error: --ncbi is needed to select the NCBI inputs to build from" >&2
    exit 1
fi

if $only_ncbi && $only_pp; then
    echo "Error: --only-ncbi and --only-pp are mutually exclusive" >&2
    exit 1
fi

if $only_ncbi; then
    run_pp_stage=false
fi

if $only_pp; then
    run_ncbi_stage=false
fi

ncbi_dir="support_files/ncbi"

# No names given: run every json under support_files/ncbi

if [ "${#ncbi_inputs[@]}" -eq 0 ]; then
    for json in "$ncbi_dir"/*.json; do
        [ -e "$json" ] || continue
        ncbi_inputs+=("$(basename "$json" .json)")
    done
fi

# Every requested name must resolve to a json in support_files/ncbi

for name in "${ncbi_inputs[@]}"; do
    if [ ! -f "$ncbi_dir/$name.json" ]; then
        echo "Error: no such NCBI input '$name' ($ncbi_dir/$name.json)" >&2
        exit 1
    fi
done

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

if $run_ncbi_stage; then
    if [ -n "$chunks" ]; then
        if ! python3 "source/ncbi/run_ncbi.py" "${ncbi_inputs[@]}" --chunks "$chunks"; then
            echo "Error: ncbi data stage failed" >&2
            exit 1
        fi
    else
        if ! python3 "source/ncbi/run_ncbi.py" "${ncbi_inputs[@]}"; then
            echo "Error: ncbi data stage failed" >&2
            exit 1
        fi
    fi
fi

if $run_pp_stage; then
    shopt -s nullglob
    pp_bins=("$data_dir"/Raw/Bins/*.fa.gz)
    shopt -u nullglob

    if [ "${#pp_bins[@]}" -eq 0 ]; then
        echo "[WARN] No bins in $data_dir/Raw/Bins; skipping the processing pipeline."
    elif ! python "source/processing_pipeline/run_pp.py"; then
        echo "Error: processing pipeline failed" >&2
        exit 1
    fi
fi

if $keep_tmps; then
    echo "Keeping tmp directory"
elif [ -d "$tmp_dir" ]; then
    if rm -rf "$tmp_dir"; then
        echo "Tmp directory removed"
    else
        echo "Error removing the tmp directory." >&2
        exit 1
    fi
fi










