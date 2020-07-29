#!/bin/bash

DIR=/storage/home/cmz5202/scratch/z.day1.0/run/ET/

find ${DIR} -name "*.nc" | sort -n > filelist.txt
