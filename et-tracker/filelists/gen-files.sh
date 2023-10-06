#!/bin/bash

# ./gen-files.sh -o files.ERA5.txt -d ~/scratch/h1files/ERA5v3/ -p "*.h1.2015*.nc"

# Defaults
OUTNAME=files.ERA5.txt
DIR=~/scratch/h1files/ERA5v3/
PATTERN="*.nc"

# Parse anything the user passed in
while getopts "o:d:p:" opt; do
    case $opt in
        o) OUTNAME="$OPTARG";;
        d) DIR="$OPTARG";;
        p) PATTERN="$OPTARG";;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

echo "DIR: $DIR"
echo "PATTERN $PATTERN"
echo "OUTNAME $OUTNAME"

find "${DIR}" -name "${PATTERN}" | sort -n > "${OUTNAME}"
