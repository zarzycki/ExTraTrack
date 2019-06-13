#!/bin/bash

TYPE="ape"
FILE=../et-tracker/climatology_files/etdetails_${TYPE}.txt

while read line; do
  echo $line
  stringarray=($line)
  itcstr=${stringarray[0]}
  yyyy=${stringarray[5]}
  mm1=`printf "%02d\n" $((10#${stringarray[6]}-10#1))`
  mm2=`printf "%02d\n" $((10#${stringarray[10]}+10#1))`
  date_start="${stringarray[5]} ${stringarray[6]} ${stringarray[7]} ${stringarray[8]}"
  date_end="${stringarray[9]} ${stringarray[10]} ${stringarray[11]} ${stringarray[12]}"
  echo ncl et_individual_nx_traj.ncl 'type="'${TYPE}'"' 'itcstr="'${itcstr}'"' 'yyyy="'${yyyy}'"' 'mm1="'${mm1}'"' 'mm2="'${mm2}'"' date_start=\""$date_start\"" date_end=\""$date_end\""
  ncl et_individual_nx_traj.ncl 'type="'${TYPE}'"' 'itcstr="'${itcstr}'"' 'yyyy="'${yyyy}'"' 'mm1="'${mm1}'"' 'mm2="'${mm2}'"' date_start=\""$date_start\"" date_end=\""$date_end\""
done < $FILE
