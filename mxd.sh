#!/bin/bash
run_layout=false
run_chen=false
run_t=false
run_all=false

if [[ "$#" -eq 1 && "$1" == "--help" ]]; then
    cat << 'EOF'
In development
EOF
    exit 0
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --layout)
            run_layout=true
            shift
            ;;
        --chen)
            run_chen=true
            shift
            ;;
        --t)
            run_t=true
            shift
            ;;
        --all)
            run_all=true
            shift
            ;;
        *)
            echo "Invalid argument: $1" >&2
            exit 1
            ;;
    esac
done

if $run_all; then
    if $run_chen || $run_t; then
        echo "Error: Argument --all was used alongside to specific subprojects, this argument enables all."
        exit 1
    fi
    run_t=true
    run_chen=true
fi

if ! $run_chen && ! $run_t; then
    echo "Error: The whole database or specific subprojects are needed to run MXD"
    exit 1
fi

echo "Starting"
echo "Creating Data Directories"

data_dir="Data"
sub_data_dirs=("$data_dir/Entities" "$data_dir/Raw" "$data_dir/Raw/Bins" "$data_dir/Raw/Processed" "$data_dir/Reports")

if [ -d "$data_dir" ]; then
    echo "Directory '$data_dir' already exists."
    read -p "Replace? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Removing existing directory..."
            if rm -rf "$data_dir"; then
                echo "Directory removed successfully."
            else
                echo "Error removing the directory."
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
    echo "Done!"
fi

if $run_chen; then
    echo "Downloading and processing the Chen et al. (2022) database"
    python3 "source/chen_data/run_chen.py"
fi

if $run_t; then
    echo "Processing the test subset"
    echo "Fingindo que rodou a padronização"
fi 

echo "Running the processing pipeline"

python3 "source/processing_pipeline/run_pp.py"








