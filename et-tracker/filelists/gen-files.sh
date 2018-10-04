#!/bin/bash

DIR=/global/homes/c/czarzyck/scratch/et-hyperion/

find ${DIR} -name "*h7*.nc" | sort -n > filelist.txt
