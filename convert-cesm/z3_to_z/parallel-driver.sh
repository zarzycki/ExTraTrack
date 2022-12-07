#!/bin/bash

#PBS -l nodes=1:ppn=16
#PBS -l walltime=12:00:00
#PBS -A open

source /storage/home/cmz5202/.bashrc

module load parallel
NUMCORES=16
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

cd /storage/home/cmz5202/sw/ExTraTrack/convert-cesm/z3_to_z

for f in $(find /storage/home/cmz5202/scratch/etc-tracking-proc/natlantic_30_x4/ -name 'natlantic30*_AMIP.cam.h6.*.nc' ); 
do
  echo $f
  echo "ncks -O -x -v Z ${f} ${f} ; ncl add-Z-to-tracker.ncl 'filename=\"'${f}'\"' ; ncks -O -4 -L 1 $f $f" >> ${COMMANDFILE}
done

parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}

