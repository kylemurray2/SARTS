#!/bin/bash

# Run this to debug an individual date when processing stack with stack processor in ISCE.
# This will make a new script to isolate all of the commands needed to do that date.
# Run the new script. 

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <date in YYYYMMDD format>"
    exit 1
fi

DATE=$1

grep "${DATE}" run_files/* | \
awk -F':' '{print $2}' | \
sed 's/&//g' > run_commands_${DATE}.sh

# Make the output script executable
chmod +x run_commands_${DATE}.sh

echo "Commands written to run_commands_${DATE}.sh"

