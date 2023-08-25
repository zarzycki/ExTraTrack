#!/bin/bash -l

################################################################
#PBS -N ET_convert_cesm
#PBS -A P93300042
#PBS -l walltime=10:49:00
#PBS -q premium
#PBS -k oe
#PBS -m a
#PBS -M zarzycki@ucar.edu
#PBS -l select=1:ncpus=36:mem=109GB
################################################################
#> #SBATCH -N 1                #Use 2 nodes
#> #SBATCH -t 07:59:00         #Set 30 minute time limit
#> #SBATCH -q regular          #Use the regular QOS
#> #SBATCH -L SCRATCH          #Job requires $SCRATCH file system
#> #SBATCH -C knl,quad,cache   #Use KNL nodes in quad cache format (default, recommended)
################################################################
#> #PBS -N interp_CESM
#> #PBS -l nodes=1:ppn=8
#> #PBS -l walltime=6:00:00
#> #PBS -A open
################################################################

starttime=$(date -u +"%s")

module load parallel

NLFILE=nl.hyperion
NUMCORES=8
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt
CASENAME=`grep -hnr "casename" ${NLFILE} | cut -d'"' -f 2`

FILES=`find /glade/u/home/zarzycki/scratch/hyperion/${CASENAME}/ -name "*cam.h2.*.nc" | sort -n`

cd /glade/u/home/zarzycki/sw/ExTraTrack/convert-cesm/se_to_extratrack

for f in ${FILES}
do
  NCLCOMMAND="ncl create-files.ncl 'f2name=\"'${f}'\"' 'nlfile=\"'${NLFILE}'\"'     "
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

# Launch GNU parallel
parallel --jobs ${NUMCORES} -u --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${COMMANDFILE}

#parallel --jobs ${NUMCORES} -u < ${COMMANDFILE}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm ${COMMANDFILE}

echo "${tottime}\n"
