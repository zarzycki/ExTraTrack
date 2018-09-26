#!/bin/bash -l

################################################################
#>#PBS -N test_gnu
#>#PBS -A P54048000 
#>#PBS -l walltime=01:49:00
#>#PBS -q premium
#>#PBS -k oe
#>#PBS -m a 
#>#PBS -M zarzycki@ucar.edu
#>#PBS -l select=1:ncpus=36:mem=109GB
################################################################
#SBATCH -N 1                #Use 2 nodes
#SBATCH -t 08:57:00         #Set 30 minute time limit
#SBATCH -q premium          #Use the regular QOS
#SBATCH -L SCRATCH          #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use KNL nodes in quad cache format (default, recommended)
################################################################

starttime=$(date -u +"%s")

module load parallel
module load ncl 

NUMCORES=12
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

STYR=1984
ENYR=2014
TYPECAPS=EXT

for DATA_YEAR in $(eval echo {$STYR..$ENYR})
do
  NCLCOMMAND="ncl reanalysis_et_cyclone_traj.NEW.ncl 'typecaps=\"${TYPECAPS}\"' 'year_min_str=\"${DATA_YEAR}\"' 'year_max_str=\"${DATA_YEAR}\"'"
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

parallel --jobs ${NUMCORES} -u < ${COMMANDFILE}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm ${COMMANDFILE}

printf "${tottime}\n" >> timing.txt


