#!/bin/bash
# Author: KM
# Date: 10/04/2023
# This script is designed for when you want to make a new time series by cropping an area of another existing stack.

# Check if an argument was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_directory>"
    exit 1
fi

# Convert relative path to absolute path
DIR_PATH=$(realpath "$1")

# Check if the given directory exists
if [ ! -d "$DIR_PATH" ]; then
    echo "Error: $DIR_PATH is not a directory!"
    exit 1
fi


# Directories to link
declare -a dirs_to_link=("baselines" "configs" "DEM" "reference" )

for dir in "${dirs_to_link[@]}"; do
    source_dir="$DIR_PATH/$dir"
    if [ -d "$source_dir" ]; then
        ln -s "$source_dir" "./$(basename "$source_dir")"
    else
        echo "Warning: $dir does not exist in $DIR_PATH!"
    fi
done

# Create a new directory named "merged" in the current directory
mkdir -p "./merged"

# Directories inside "merged" to link
declare -a dirs_to_link=("SLC" "baselines")

for dir in "${dirs_to_link[@]}"; do
    source_dir="$DIR_PATH/merged/$dir"
    if [ -d "$source_dir" ]; then
        ln -s "$source_dir" "./merged/"
    else
        echo "Warning: $dir does not exist in $DIR_PATH/merged/"
    fi
done

# Handle geom_reference directory in "merged"
geom_ref_source="$DIR_PATH/merged/geom_reference"
geom_ref_dest="./merged/geom_reference"

# Create the geom_reference directory inside the "merged" directory in the current directory
mkdir -p "$geom_ref_dest"

# Link files from geom_reference that don't contain 'crop' or 'lk' in their titles
for file in "$geom_ref_source"/*; do
    file_basename=$(basename "$file")
    if [[ ! "$file_basename" =~ (crop|lk) ]]; then
        ln -s "$file" "$geom_ref_dest/$file_basename"
    fi
done

# Copy the file "params.yaml" from the original given directory to the current directory
params_file="$DIR_PATH/params.yaml"
if [ -e "$params_file" ]; then
    cp "$params_file" .
    echo "params.yaml has been copied to the current directory."
else
    echo "Error: params.yaml does not exist in $DIR_PATH!"
    exit 1
fi
