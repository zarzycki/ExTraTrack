#!/bin/bash

# Driver script to search for ETCs in unstructured CAM-SE data
# Colin Zarzycki, 3/21/16

############ USER OPTIONS #####################

## Path to TempestExtremes binaries on YS
TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes_noMPI/

UQSTR=MP15
TOPOFILE=/glade/u/home/zarzycki/scratch/TEST-ET3/f.asd2017.cesm20b05.FAMIPC6CLM5.mp15a-120a-US_t12.exp213.cam.TOPO.nc
PATHTOFILES=/glade/scratch/zarzycki/TEST-ET3/
CONNECTFLAG="" 

FILES=`ls ${PATHTOFILES}/*.cam.h7.*.nc`

############ TRACKER MECHANICS #####################

# Loop over files to find candidate cyclones
#rm cyc.${UQSTR} trajectories.txt.${UQSTR}
#touch cyc.${UQSTR}
STR_DETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyc_tempest.${UQSTR} --closedcontourcmd PSL,300.0,4.0,0;_DIFF(Z(0),Z(5)),-10.0,5,0.25 --mergedist 5.0 --searchbymin PSL --outputcmd PSL,min,0;_VECMAG(UBOT,VBOT),max,2;PHIS,max,0"
for f in ${FILES[@]};
do
  echo "Processing $f..."
  ${TEMPESTEXTREMESDIR}/bin/DetectCyclonesUnstructured --in_data "${f};${TOPOFILE}" ${STR_DETECT}
  cat cyc_tempest.${UQSTR} >> cyc.${UQSTR}
  rm cyc_tempest.${UQSTR}
done

# Stitch candidate cyclones together
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "j,i,lon,lat,slp,wind,phis" --range 5.0 --minlength 8 --maxgap 2 --in cyc.${UQSTR} --out trajectories.txt.${UQSTR} --threshold "phis,<=,500,8;lat,<=,45.0,8;lat,>=,-45.0,8"

#rm cyc.${UQSTR}   #Delete candidate cyclone file
