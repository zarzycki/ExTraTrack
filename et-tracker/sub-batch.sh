#!/bin/bash -l

#PBS -N ExTraTrack
#PBS -A UPSU0032
#PBS -l walltime=01:49:00
#PBS -q premium
#PBS -k oe
#PBS -m a
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB

#>#SBATCH -N 1                #Use 2 nodes
#>#SBATCH -t 09:57:00         #Set 30 minute time limit
#>#SBATCH -q premium          #Use the regular QOS
#>#SBATCH -L SCRATCH          #Job requires $SCRATCH file system
#>#SBATCH -C haswell   #Use KNL nodes in quad cache format (default, recommended)

starttime=$(date -u +"%s")

DOPOSTPROCESS=true

module load parallel
module load ncl

NUMCORES=18
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

STYR=1985
ENYR=2014
NAMELISTFILE="./user-nl/nl.hyp"

for DATA_YEAR in $(eval echo {$STYR..$ENYR})
do
  NCLCOMMAND="ncl ExTraTrack.ncl 'nlfile=\"${NAMELISTFILE}\"' 'year_min_str=\"${DATA_YEAR}\"' 'year_max_str=\"${DATA_YEAR}\"'"
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

if [[ $HOST = *"cori"* ]]; then
  parallel --jobs ${NUMCORES} -u < ${COMMANDFILE}
else
  parallel --jobs ${NUMCORES} -u --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${COMMANDFILE}
fi

if [ "$DOPOSTPROCESS" = true ] ; then
  ncl et_concat_trajs.ncl 'nlfile="'${NAMELISTFILE}'"'
  ncl et_avg_text.ncl 'nlfile="'${NAMELISTFILE}'"'
  ncl et_yearly_clim.ncl 'nlfile="'${NAMELISTFILE}'"'
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm ${COMMANDFILE}

printf "${tottime}\n" >> timing.txt


