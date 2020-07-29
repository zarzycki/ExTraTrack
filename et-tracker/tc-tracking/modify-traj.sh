#!/bin/bash
  
NEWFILE=tc_traj
rm $NEWFILE

STORMID=0
while IFS='' read -r line || [[ -n "$line" ]]; do
#  echo $line
  if [ "${line:0:1}" == "s" ] ; then
    echo $line
    STORMSTRING=`printf "%03d\n" $STORMID`
    echo ${line}" "${STORMSTRING} >> ${NEWFILE}
    STORMID=$((STORMID+1))
  else
    eval x=($line)

    #parse MSLP
    var=`echo ${x[4]} | sed -e 's/[eE]+*/\\*10\\^/'`
    mslp=$(bc -l <<< "scale=3;${var}/100")
    #parse wind
    var=`echo ${x[5]} | sed -e 's/[eE]+*/\\*10\\^/'`
    wind=$(bc -l <<< "scale=3;${var}/1")

#        314     145     157.000000      -17.500000      9.964526e+04    1.532690e+01    0.000000e+00    1999    1       20      12

    # These lines have to be modified based on how your columns work out (i.e., unstruct, do you have PHIS, etc.)

    #echo ${x[2]}' '${x[3]}' '${wind}' '${mslp}' '${x[7]}' '${x[8]}' '${x[9]}' '${x[10]}  >> ${NEWFILE}
    echo ${x[2]}' '${x[3]}' '${wind}' '${mslp}' '${x[6]}' '${x[7]}' '${x[8]}' '${x[9]}  >> ${NEWFILE}

  fi
done < "$1"
