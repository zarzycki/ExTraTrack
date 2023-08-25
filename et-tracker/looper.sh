#!/bin/bash

ARR_STYR=("1985" "1985" "1985" "2070" "2070" "2070" "2070" "2070" "2070")
ARR_ENYR=("2014" "2014" "2014" "2100" "2100" "2100" "2100" "2100" "2100")
ARR_NLS=("nl.hyp.ref.dtime900" "nl.hyp.ref.dtime900.002" "nl.hyp.ref.dtime900.003" "nl.hyp.ref.dtime900.rcp45" "nl.hyp.ref.dtime900.rcp45.002" "nl.hyp.ref.dtime900.rcp45.003" "nl.hyp.ref.dtime900.rcp85" "nl.hyp.ref.dtime900.rcp85.003" "nl.hyp.ref.dtime900.rcp85.004")

length=${#ARR_NLS[@]}

for ((i=0; i<$length; i++)); do
  echo "Loop: $i"
  sed -i "s?^STYR=.*?STYR=${ARR_STYR[$i]}?" sub-batch.sh
  sed -i "s?^ENYR=.*?ENYR=${ARR_ENYR[$i]}?" sub-batch.sh
  sed -i "s?^NAMELISTFILE=.*?NAMELISTFILE=\"./user-nl/${ARR_NLS[$i]}\"?"  sub-batch.sh
  qsub sub-batch.sh
  sleep 5
done