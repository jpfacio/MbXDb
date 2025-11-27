#!/bin/bash

bins="Data/Chen_Analysis/Bins"
unzip_dir="Data/Chen_Analysis/Bins/TMP_UNZIP"
output="Data/Chen_Analysis/CheckM2"

if [ ! -d "$bins" ]; then
    echo "Error: Directory '$bins' not found."
    exit 1
fi

mkdir -p "$unzip_dir"

echo "Creating temp diretory for checkm2 analysis"

for i in "$bins"/*.fa.gz; do
    if [ -f "$i" ]; then
        gzip -dc "$i" > "$unzip_dir/$(basename "$i" .gz)"
    fi
done

if ! ls "$unzip_dir"/*fa >/dev/null 2>&1; then
    echo "Unzipping error"
    rm -rf "$unzip_dir"
    exit 1
fi

echo "Running checkm2.."

mkdir -p "$output"

checkm2 predict \
    --input "$unzip_dir" \
    --output-directory "$output" \
    --threads 8 \
    -x fa \
    --force

rm -rf "$unzip_dir"

echo "Concluded!"


