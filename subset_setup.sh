#!/bin/bash

# Author: KM
# Date: 10/04/2023

# This script is designed for when you want to make a new time series by cropping an area of another existing stack.
#  You make a new working directory.  Then run this script in that directory.
#  It will link to the key directories from the main stack, and copy the parmas.yaml file.
#  Only the full geom files are linked so that new cropped and downlooked geom files can be added without overwriting the main stack files.

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

# Create a new directory named "merged" in the current directory
mkdir -p "./merged"

# Directories to link
declare -a dirs_to_link=("baselines" "DEM" "reference" )

for dir in "${dirs_to_link[@]}"; do
    source_dir="$DIR_PATH/$dir"
    if [ -d "$source_dir" ]; then
        ln -s "$source_dir" "./$(basename "$source_dir")"
    else
        echo "Warning: $dir does not exist in $DIR_PATH!"
    fi
done

# Directories inside "merged" to link
declare -a dirs_to_link=( "baselines")

for dir in "${dirs_to_link[@]}"; do
    source_dir="$DIR_PATH/merged/$dir"
    target_link="./merged/$dir"
    
    if [ -d "$source_dir" ]; then
        if [ ! -e "$target_link" ]; then
            ln -s "$source_dir" "./merged/"
        else
            echo "Warning: $target_link already exists. Skipping linking."
        fi
    else
        echo "Warning: $dir does not exist in $DIR_PATH/merged/"
    fi
done

# Define the source and target base directories
src_dir="$DIR_PATH/merged"
target_dir="./merged/SLC"

# Create the target base directory if it doesn't exist
mkdir -p "$target_dir"

# Iterate over each subdirectory in the source directory
for subdir in "$src_dir"/*; do
    if [ -d "$subdir" ]; then
        # Extract the name of the subdirectory
        subdir_name=$(basename "$subdir")

        # Create a corresponding subdirectory in the target directory
        mkdir -p "$target_dir/$subdir_name"

        # Link all files that do not contain 'crop' in their names
        for file in "$subdir"/*; do
            if [[ ! $(basename "$file") == *crop* ]]; then
                ln -s "$PWD/$file" "$target_dir/$subdir_name"
            fi
        done
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
    if [[ ! "$file_basename" =~ (crop|lk|water|land) ]]; then
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

echo ' '
echo 'Edit crop bounds of params.yaml'
echo 'Change Land cover file, if using water mask'