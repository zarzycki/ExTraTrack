# ExTraTrack v0.1

## aka the Extratropical Transition (tropical cyclone) Tracker

Colin M. Zarzycki and Diana R. Thatcher.

NOTE: This is pre-release code. It mainly consists of NCL / Bash scripting. It was used to track the extratropical transition of tropical cyclones in high-resolution CESM data as well as various reanalysis productions. Porting to a more flexible language (Python?) and generalizing the code is on my to-do list, but since I am no longer funded on this project, pretty far down the list. This may be an area of future work, but as of now the code is not being actively developed (save for a few cleanings, etc. when necessitated by collaborators using the code). Please e-mail me if you are interested in adapting this code for use with your products!

WARNING: This code has not been extensively verified beyond the particular projects noted below. If you find any bugs or inconsistencies, please contact zarzycki@ucar.edu.

This code is used in the paper:
Zarzycki, C. M., D. R. Thatcher, and C. Jablonowski (2017), Objective tropical cyclone extratropical transition detection in high-resolution reanalysis and climate model data, J. Adv. Model. Earth Syst., 9, 130–148, doi:10.1002/2016MS000775.

It applies cyclone phase space as defined by:
Hart, R.E., 2003: A Cyclone Phase Space Derived from Thermal Wind and Thermal Asymmetry. Mon. Wea. Rev., 131, 585–616, doi:10.1175/1520-0493(2003)131<0585:ACPSDF>2.0.CO;2 

To use this code, you *must* have gridded reanalysis or climate model data that contains variables described below.

**A sample of gridded data + TC trajectories to start with can be downloaded from:** http://www.cgd.ucar.edu/staff/zarzycki/files/ERA-sample-ETC-tracker-data.gz

The code reads in TC trajectories defined by a TC tracker, such as that defined in:
* Zarzycki and Jablonowski, JAMES, 2014
* Ullrich and Zarzycki, GMD, 2017
* Zarzycki and Ullrich, GRL, 2017
* ... among many others.

Given those TC trajectories, it searches gridded data matching the period of analysis. During the times a TC trajectory exists, it calculates Vut, Vlt, and B from Hart et al., (2003).

Following the termination of a defined warm-core TC trajectory, the code continues following features based on sea level pressure minima at subsequent time slices until the cyclone dissipates or leaves the domain.

## General procedure

1. Generate TC tracks using TC software
2. Post-process TC tracks into correct format
3. Modify `reanalysis_et_cyclone_traj.ncl` and/or gridded netCDF files to allow for variables to be read by `reanalysis_et_cyclone_traj.ncl` correctly.
4. Run `reanalysis_et_cyclone_traj.ncl` to extract extended tracks with CPS variables (B, VUT, VLT)
  * Concatenate ET trajectories output from `reanalysis_et_cyclone_traj.ncl`
5. (optional, but highly recommended) run `et_avg.ncl` to "smooth" CPS parameters with time.
6. Process/plot ET climatological statistics.

## Detailed procedure

*** Note, actual ET tracker is in ./et-tracker.
*** Functions associated with ET tracker are in ./functions
*** Plotting scripts (if available) are in ./plotting-scripts

A TC tracker (such as TempestExtremes or that used by GFDL) needs to be used to generate/initialize TC tracks. This file will be shorthanded `${TCTRAJ}` in this README.

Tracks *must* be post-processed to be in this text format

```
start 12 1980 5 1 0 0000
260.156342 6.666654 10.399 1005.583 1980 5 1 0
260.156342 7.368407 12.425 1007.042 1980 5 1 6
259.453217 8.070160 12.061 1005.444 1980 5 1 12
259.453217 8.771913 8.928 1008.520 1980 5 1 18
258.750092 8.771913 9.844 1005.786 1980 5 2 0
258.750092 9.473666 8.641 1007.734 1980 5 2 6
258.046967 9.473666 8.332 1005.352 1980 5 2 12
258.046967 9.473666 9.388 1007.724 1980 5 2 18
258.046967 9.473666 10.366 1005.119 1980 5 3 0
257.343842 10.175419 10.204 1006.493 1980 5 3 6
257.343842 10.175419 11.325 1003.938 1980 5 3 12
257.343842 10.175419 11.262 1007.239 1980 5 3 18
start 32 1980 8 3 18 0001
302.343842 12.982431 14.584 1006.678 1980 8 3 18
300.937592 12.982431 11.796 1006.424 1980 8 4 0
299.531342 12.982431 12.673 1005.519 1980 8 4 6
297.421967 14.385937 15.288 1006.124 1980 8 4 12
...
```

each trajectory consists of a header line...
```
start NUMOBS STYYYY STMM STDD STHH STORMID
```

each line of the trajectory (there should be NOBS of them) are...
```
LON LAT BOT_WINDMAG PSL YYYY MM DD HH
```

where...
```
PSL          ; minimum pressure value (Pa)
BOT_WINDMAG  ; maximum wind speed at 10m or in lowest model level (m/s)
```



`${TCTRAJ}` is usually placed in `./text_files/` and is used by `reanalysis_et_cyclone_traj.ncl` and `et_yearly_clim.ncl`

The individual trajectories with B, Vut, and Vlt can be calculated by invoking
```$> ncl reanalysis_et_cyclone_traj.ncl 'year_min_str="'${YYYY}'"' 'year_max_str="'${YYYY}'"'```

where `reanalysis_et_cyclone_traj.ncl` points to `${TCTRAJ}`. Other user settings are described in the script. `year_min` and `year_max` can be used to parallelize the code (i.e., for a 20 year dataset) you could spawn 20 single-core jobs. Most helpful for high-resolution data with many TCs.

| Namelist Key | Namelist sample | Description |
| --- | --- | --- |
| type | type="era", | Shortcode for file naming |
| tfile | tfile="./tc-tracking/tc\_traj", | Full/relative path to `${TCTRAJ}` |
| filelist | filelist="./filelist.txt", | Path to gridded netCDF files, one per line and chronologically sorted |
| basin | basin=1, | Only calculate CPS values for specific basin (< 0 = all basins) |
| latmin | latmin=-80.0, | Minimum latitude to extract from gridded data |
| latmax | latmax=80.0, | Maximum latitude to extract from gridded data |
| lonmin | lonmin=-180.0, | Minimum longitude to extract from gridded data<sup>[1](#namelistfoot1)</sup> |
| lonmax | lonmax=360.0, | Maximum latitude to extract from gridded data<sup>[1](#namelistfoot1)</sup> |
| hrintvl | hrintvl=6.0, | Time interval (hours) between data points |
<a name="namelistfoot1">1</a>: If using global data, -180. -> 360. allows for all lon orderings. For more regional definitions, this will be specific to the gridded dataset coordinate conventions.

REQUIRED DATA: Currently, variables/names needed by `reanalysis_et_cyclone_traj.ncl` are:
1. PSL (sea level pressure)
2. UBOT (lowest model level U wind)
3. VBOT (lowest model level V wind)
4. Z (geopotential height) @ 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 
    875, 900, 925, 950, 975, 1000 mb.

PSL, UBOT, and VBOT are 2-D variables at each time.
Z is 3-D, with the vertical dimension named nlev.

Example ERA-I files are included in `./sample_data/`

See block of text in `reanalysis_et_cyclone_traj.ncl` to modify data ingestion for different formats.

`reanalysis_et_cyclone_traj.ncl` will output new, storm-by-storm trajectory files in `./text_files/`

Each file will look like...
```
start   051  1980      08    03    18    001
   -57.66   12.98   1006.68   14.6   -999.00   -999.00   -999.00   -999.00   -999.00   1980   8   3  18
   -59.06   12.98   1006.42   11.8    152.38    270.00     -2.41    -11.16     44.55   1980   8   4   0
   -60.47   12.98   1005.52   12.7    152.38    270.00     -2.51     21.88     45.18   1980   8   4   6
...
```

header format:
```
start NUMOBS STYYYY STMM STDD STHH STORMID
```

each line (NOBS of them) will be...
```
LON LAT PSL BOT_WINDMAG DIST ANG B VLT VUT YYYY MM DD HH
```

where...

```
PSL          ; minimum pressure value (Pa)
BOT_WINDMAG  ; maximum wind speed at 10m or in lowest model level (m/s)
DIST         ; distance traveled between each time step (m)
ANG          ; angle of storm travel (deg)
B            ; B parameter (Hart 2003)
VLT          ; lower troposphere thermal wind (Hart 2003)
VUT          ; upper troposphere thermal wind (Hart 2003)
```


In `./text_files/` the individual `tmp_XXX_NNN.txt` files need to be concatenated into a singlular traj file
`$> cat tmp_era_* > ${ETCTRAJ_ORIG}`

You can then run 
```$> ncl et_avg_text.ncl```
... to produce a "smoothed" ET trajectory file from `./text_files/${ETCTRAJ_ORIG}`.
Settings are contained in the user options at the top of the script.
This produces a file `./text_files/${ETCTRAJ_AVG}`.

Statistics can be calculated by running:
```$> ncl et_yearly_clim.ncl```

This will produce a series of files in `./climatology_files/`
Here we assume ERA-interim 1980-2002 (as in Zarzycki et al., 2017).

*********************************************
`monthly_era_1980_2002.txt`

```
1980  01   00   00   00   00   00   00
```

```
YYYY  MM  MATL MTC  MWC  MCC  MET  MNO
```

* MATL number of atlantic basin storms
* MTC  number of storms that remained tropical
* MWC  number of storms that dissipated as warm cores
* MCC  number of storms that dissipated as cold cores
* MET  number of storms that complete transition
* MNO  number of storms that complete partial transition

*********************************************
`yearly_era_1980_2002.txt`

```
1980   02   02   00   00   00   00
```

```
YYYY   MATL MTC  MWC  MCC  MET  MNO
````

* MATL number of atlantic basin storms
* MTC  number of storms that remained tropical
* MWC  number of storms that dissipated as warm cores
* MCC  number of storms that dissipated as cold cores
* MET  number of storms that complete transition
* MNO  number of storms that complete partial transition

*********************************************
`etpath_yearly_era_1980_2002.txt`

```
1980   00   00   00
```

```
YYYY  PATH1 PATH2 PATH3
```

* PATH1 - storms undergoing Type I ET
* PATH2 - storms undergoing Type II ET
* PATH3 - storms undergoing Type III ET
* (see Zarzycki et al., 2017 for type definitions)

*********************************************
`storms_era.txt`

(list of storms, 1 per line, that undergo total ET)

*********************************************
`life_era.txt`

```
storm    0005   028    1981  08  19  00
         0005   005    1981  08  20  12    1981  08  21  18
```

```
storm   STORMID   TOTDUR  TCSTARTYYYYMMDDHH
         STORMID ETDUR ETSTARTYYYMMMDDHH  ETENDYYYMMMDDHH
```

* TOTDUR - total length of TC + ETC tracked
* ETDUR - time between ET onset and ET completion
* TCSTARTYYYYMMDDHH - initial time of TC tracking
* ETSTARTYYYMMMDDHH - time of ET onset
* ETENDYYYMMMDDHH - time of ET completion

Note: TOTDUR/ETDUR in "timesteps" (number of datapoints, so 5 for 6-hourly pts = 30 hours)

*********************************************
`etdetails_era.txt`

```
 0005   005    1  1002.36  999.80   1981  08  20  12    1981  08  21  18
```

```
STORMID ETDUR PATH SLPST   SLPEN     ETSTARTYYYMMMDDHH   ETENDYYYMMMDDHH
```

* PATH is path type defined in Zarzycki et al., 2017 
* SLPST - SLP at start of transition
* SLPEN - SLP at end of transition
