#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 1                      # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J sub_ncl 
#BSUB -W 3:59                  # wall clock limit
#BSUB -P P05010048               # account number

################################################################

date

module swap ncl ncl/6.3.0

TYPECAPS=ERA
YYYY=2001
ncl reanalysis_et_cyclone_traj.ncl 'typecaps="'${TYPECAPS}'"' 'year_min_str="'${YYYY}'"' 'year_max_str="'${YYYY}'"'

date 
