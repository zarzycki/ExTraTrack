#!/bin/bash -l

#PBS -N ExTraTrack
#PBS -A P93300642
#PBS -l walltime=23:30:00
#PBS -q casper@casper-pbs
#PBS -k oe
#PBS -m a
#PBS -l select=1:ncpus=10:mem=500GB

#>#SBATCH -N 1                #Use 2 nodes
#>#SBATCH -t 09:57:00         #Set 30 minute time limit
#>#SBATCH -q premium          #Use the regular QOS
#>#SBATCH -L SCRATCH          #Job requires $SCRATCH file system
#>#SBATCH -C haswell   #Use KNL nodes in quad cache format (default, recommended)

starttime=$(date -u +"%s")

DOPOSTPROCESS=true

module load parallel
module load ncl
module load peak_memusage

NUMCORES=10
TIMESTAMP=`date +%s%N`
TIMESTAMP="1696429188031086934"
COMMANDFILE=commands.${TIMESTAMP}.txt
YYYYMMDDHHSS=`date '+%Y%m%d%H%M'`

STYR=1980
ENYR=2018
NAMELISTFILE="./user-nl/nl.ERA5"

for DATA_YEAR in $(eval echo {$STYR..$ENYR})
do
  NCLCOMMAND="ncl ExTraTrack.ncl 'nlfile=\"${NAMELISTFILE}\"' 'year_min_str=\"${DATA_YEAR}\"' 'year_max_str=\"${DATA_YEAR}\"' 'UQSTR=\"${TIMESTAMP}\"'  "
  echo "${NCLCOMMAND} >> logs/log.${YYYYMMDDHHSS}.${PBS_JOBID}.${DATA_YEAR}" >> ${COMMANDFILE}
done

mkdir -p logs/

# build lookup table
#ncl et_build_lookup.ncl 'nlfile="'${NAMELISTFILE}'"' 'UQSTR="'${TIMESTAMP}'"'

if [[ $HOST = *"cori"* ]]; then
  parallel --jobs ${NUMCORES} -u < ${COMMANDFILE}
else
  peak_memusage.exe parallel --jobs ${NUMCORES} --workdir $PWD < ${COMMANDFILE}
fi

if [ "$DOPOSTPROCESS" = true ] ; then
  ncl et_concat_trajs.ncl 'nlfile="'${NAMELISTFILE}'"'
  ncl et_avg_text.ncl 'nlfile="'${NAMELISTFILE}'"'
  ncl et_yearly_clim.ncl 'nlfile="'${NAMELISTFILE}'"'
  ncl et_to_te_file.ncl 'nlfile="'${NAMELISTFILE}'"'
fi

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm -fv ${COMMANDFILE}
rm -fv lookup_${TIMESTAMP}.nc

printf "${tottime}\n" >> timing.txt


