#!/bin/bash

### 8/29/17 - zarzycki@ucar.edu
### This script concatenates the traj text files generated by reanalysis_et_cyclone_traj.ncl
### It is obviously not very complicated :)

type="era"
cd ./text_files
cat tmp_${type}_*.txt > traj_et_${type}_orig
#rm tmp_${type}_*.txt

