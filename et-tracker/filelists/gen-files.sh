#!/bin/bash

OUTNAME=/glade/u/home/zarzycki/sw/ExTraTrack/et-tracker/filelists/files.dtime900.002.txt
DIR=/glade/u/home/zarzycki/scratch/hyperion/ET2_CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.002

find ${DIR} -name "*.nc" | sort -n > ${OUTNAME}
