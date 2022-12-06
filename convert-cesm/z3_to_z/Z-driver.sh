#!/bin/bash

# Given a CAM NetCDF file with Z3 in hybrid coordinates, create a p_lev version of Z
# In other words, add Z on constant pressure surfaces which is suitable for ExTraTrack

ARCHIVEDIR=~/scratch/etc-tracking-proc/natlantic_30_x4/1999/
FILES=`ls ${ARCHIVEDIR}/*.h6.*.nc`

for f in $FILES
do
  echo $f
  ncks -O -x -v Z ${f} ${f}
  ncl add-Z-to-tracker.ncl 'filename="'${f}'"'
done
