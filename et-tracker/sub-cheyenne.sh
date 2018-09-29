#!/bin/bash -l

################################################################
#PBS -N test_gnu
#PBS -A P54048000 
#PBS -l walltime=01:49:00
#PBS -q premium
#PBS -k oe
#PBS -m a 
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB
################################################################

starttime=$(date -u +"%s")

module load parallel
module load ncl 

NUMCORES=3
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

STYR=2010
ENYR=2012
TYPECAPS=EXT

for DATA_YEAR in $(eval echo {$STYR..$ENYR})
do
  NCLCOMMAND="ncl reanalysis_et_cyclone_traj.ncl 'typecaps=\"${TYPECAPS}\"' 'year_min_str=\"${DATA_YEAR}\"' 'year_max_str=\"${DATA_YEAR}\"'"
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

parallel --jobs ${NUMCORES} -u --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${COMMANDFILE}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm ${COMMANDFILE}

printf "${tottime}\n" >> timing.txt


