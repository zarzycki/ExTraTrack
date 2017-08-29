#!/bin/bash

## This script is designed to spawn NUMYEARS instances of some script
## which is initiated for a given year

STYR=1980
ENYR=2002
TYPECAPS=ERA

thisdir=$PWD
scriptdir=${thisdir}
scriptname=sub-geyser.sh

for DATA_YEAR in $(eval echo {$STYR..$ENYR})
do
  echo "Sedding ${DATA_YEAR}"
  sed -i "s?^YYYY.*?YYYY=${DATA_YEAR}?" ${scriptdir}/${scriptname}
  sed -i "s?^TYPECAPS.*?TYPECAPS=${TYPECAPS}?" ${scriptdir}/${scriptname}

  cd ${scriptdir}
  LSB_JOB_REPORT_MAIL=N bsub < ${scriptname}
  cd $thisdir
  sleep 3
    
done
