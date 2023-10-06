#!/bin/bash

ARR_CONFIG=("dtime900" "dtime900.002" "dtime900.003" "dtime900.rcp45" "dtime900.rcp45.002" "dtime900.rcp45.003" "dtime900.rcp85" "dtime900.rcp85.003" "dtime900.rcp85.004")
length=${#ARR_CONFIG[@]}

for ((i=0; i<$length; i++)); do
  CONFIG=${ARR_CONFIG[$i]}
  ncl et_to_te_file.ncl 'nlfile="user-nl/nl.hyp.ref.'${CONFIG}'"'
done
