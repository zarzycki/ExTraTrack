#!/bin/bash -l

################################################################
#PBS -N test_gnu
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

NUMCORES=8
TIMESTAMP=`date +%s%N`
COMMANDFILE=commands.${TIMESTAMP}.txt

FILES=`find /glade/u/home/zarzycki/scratch/hyperion/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900/atm/ -name "*cam.h2.*.nc" | sort -n`

cd /glade/u/home/zarzycki/sw/ExTraTrack/convert-cesm/se_to_extratrack

for f in ${FILES}
do
  NCLCOMMAND="ncl create-files.ncl 'f2name=\"'${f}'\"' 'nlfile=\"nl.hyperion\"'     "
  echo ${NCLCOMMAND} >> ${COMMANDFILE}
done

# Launch GNU parallel
parallel --jobs ${NUMCORES} -u --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < ${COMMANDFILE}

#parallel --jobs ${NUMCORES} -u < ${COMMANDFILE}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))

rm ${COMMANDFILE}

echo "${tottime}\n"


