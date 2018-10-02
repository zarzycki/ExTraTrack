#!/bin/bash

DIR=/global/homes/c/czarzyck/scratch/et-hyperion/

find ${DIR} -name "*h7.1984*.nc" | sort -n > filelist.txt
