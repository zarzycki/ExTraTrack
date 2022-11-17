# ExTraTrack

## aka the Extratropical Transition (tropical cyclone) Tracker

Colin M. Zarzycki and Diana R. Thatcher.

ExTraTrack allows for the calculation of cyclone phase space (CPS) parameters in gridded reanalysis/climate data given a set of cyclone trajectories. The code reads in trajectory data from standard tropical cyclone trackers (or observed trajectories for reanalysis) and calculates thermal symmetry and warm/cold core depth from the geopotential field by collocating the cyclone center (defined by PSL min) in space and time. The tracker can follow PSL minima following the termination of TC trajectories to fully encompass storm lifecycles.

This code is used in the paper:
_Zarzycki, C. M., D. R. Thatcher, and C. Jablonowski (2017), Objective tropical cyclone extratropical transition detection in high-resolution reanalysis and climate model data, J. Adv. Model. Earth Syst., 9, 130–148, doi:10.1002/2016MS000775._

**NOTE:** This is pre-release code that mainly consists of NCL / Bash scripting. Porting to a more flexible language (Python?) and generalizing the code is on my to-do list, but not a high priority. Please e-mail me if you are interested in adapting this code for use with your products. If you are a user and find a bug or think you can contribute an improvement, open a ticket, pull request, etc.

**WARNING:** This code has not been extensively verified beyond the particular projects noted below. If you find any bugs or inconsistencies, please contact czarzycki@psu.edu.

It applies cyclone phase space as defined by:
_Hart, R.E., 2003: A Cyclone Phase Space Derived from Thermal Wind and Thermal Asymmetry. Mon. Wea. Rev., 131, 585–616_

To use this code, you *must* have gridded reanalysis or climate model data that contains variables described below.

**A sample of gridded data + TC trajectories to start with can be downloaded from:** [http://www.colinzarzycki.com/files/ERA-sample-ETC-tracker-data.tar.gz](http://www.colinzarzycki.com/files/ERA-sample-ETC-tracker-data.tar.gz)

## General procedure

1. Generate TC trajectories using TC software
2. 'Massage' gridded netCDF files to standardized, CF-compliant, ExTraTrack-supported format
3. Build "filelist" from netCDF files
4. Run `ExTraTrack.ncl` to extract extended tracks with CPS variables (B, VUT, VLT)
5. Concatenate ET trajectories output, "smooth" CPS parameters, process/plot ET climatological statistics

## Detailed procedure

*** Note, actual ET tracker is in ./et-tracker.
*** Functions associated with ET tracker are in ./functions
*** Plotting scripts (if available) are in ./plotting-scripts

### 1.) Generate TC tracks
A TC tracker (such as [TempestExtremes](https://github.com/ClimateGlobalChange/tempestextremes) or TECA) needs to be used to generate/initialize TC tracks. This file will be shorthanded `${TCTRAJ}` in this README.

For more information about TempestExtremes, please see...
* Ullrich and Zarzycki, GMD, 2017
* Zarzycki and Ullrich, GRL, 2017

Sample scripts for generating TC trajectories using TempestExtremes are found in `./tc-tracking/`

`${TCTRAJ}` is usually placed in `./text_files/` and is used by `reanalysis_et_cyclone_traj.ncl` and `et_yearly_clim.ncl`

**NOTE:** ExTraTrack by default wraps over the top of the TempestExtremes track format. Files are ASCII containing multiple trajectories. They can be delimited by tabs or spaces.

each trajectory consists of a header line...
```
start NUMOBS STYYYY STMM STDD STHH
```

each line of the trajectory (there should be NOBS of them) are...
```
IX JY LON LAT WIND PRES VAR1 VAR2 ... VARN PSL YYYY MM DD HH
```

**For ExTraTrack to parse this file correctly, the user must tell ExTraTrack what columns (0-based) correspond to the lon, lat, maximum surface wind, and minimum pressure values. This is specified by namelist parameter:

```
trajinds=2,3,4,5,
```

For example, the above tells ExTraTrack that the lon is in the third column, lat is in the fourth, etc. (again, zero indexed).

More information is below.

### 2.) Produce gridded, regular lat-lon netCDF

Variables must be in regular lat-lon format with dimensions ntime x nlev x nlat x nlon.

REQUIRED DATA: Currently, variables/names needed by `reanalysis_et_cyclone_traj.ncl` are:
1. PSL (sea level pressure in Pa)
2. UBOT (lowest model level U wind, in m/s)
3. VBOT (lowest model level V wind, in m/s)
4. Z (geopotential height, in m) @ 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850,
    875, 900, 925, 950, 975, 1000 mb.

PSL, UBOT, and VBOT are 2-D variables at each time. (ntime x nlat x nlon)
Z is 3-D, with the vertical dimension named "lev." (ntime x nlev x nlat x nlon)

"lat" must be increasing from south (-90) to north (90) (degrees\_north) and may be a subset of the global domain
"lon" must be increasing from west to east (degrees\_east). The preferred convention is 0-360, but -180-180 should be correctly handled as well.
"lev" must be in mb (hPa), from top-to-bottom, and lev must be an associated coordinate of the Z array.

The time dimension must be the record dimension and have a CF-compliant units attribute (i.e., "days since 2000-01-01"). There should also be a calendar attribute if using a calendar other than "noleap"

Files may only be split along the time dimension. In other words, PSL, UBOT, VBOT, and Z at time t must be on the same file, but files can be split to have a single time per file, 4 times per file, 1460 times per file, etc. Any files with more than one month of data will occur additional memory overhead, so daily/weekly/monthly are preferred.

Example ERA-I files that are compliant are included in the example tar.gz file listed above.

Users may modify the data ingestion within `reanalysis_et_cyclone_traj.ncl` but only the above format will be supported.

### 3.) Generate a file list

A sorted (from oldest data date to newest), one-per-line list of absolute paths to the processed netCDF files must be generated. The simplest way to do so is to use UNIX's `find` utility. For example, if I have all of my files stored under `/home/$LOGNAME/mydata/` in an arbitrary directory structure, I can use...

`find /home/$LOGNAME/mydata/ -name "*.nc" | sort -n > myfilelist.txt`

This filelist will need to be referenced in the namelist described below.

### 4.) Calculate CPS with NCL code

The individual trajectories with B, Vut, and Vlt can be calculated by invoking

```$> ncl ExTraTrack.ncl```

**NOTE:** The code requires reading a namelist, the default option is to read the file titled `namelist` from the same directory as ExTraTrack.

`namelist` defines an external namelist with specific configuration options. A sample namelist is included in the example tarball. Keys must be defined as `VAR="TEXT",` when a string or `VAR=999.9` when numeric.

| Namelist Variable | Namelist sample | Type | Description |
| --- | --- | --- | --- |
| type | type="era", | string | Shortcode describing particular dataset |
| tfile | tfile="./tc-tracking/tc\_traj", | string | Path to `${TCTRAJ}`. Can be relative, but safer as full path. |
| filelist | filelist="./filelist.txt", | string | Path to gridded netCDF files, one per line and chronologically sorted (see step 4)|
| etfileori | etfileori="./text\_files/traj\_et\_era\_orig", | string | Concatenated ET file (no smoothing) |
| etfileavg | etfileavg="./text\_files/traj\_et\_era\_avg", | string | Concatenated ET file (post smoothing)|
| basin | basin=1, | int | Only calculate CPS values for specific basin (set to <0 for all basins) |
| latmin | latmin=-80.0, | float | Minimum latitude to extract from gridded data |
| latmax | latmax=80.0, | float | Maximum latitude to extract from gridded data |
| lonmin | lonmin=-180.0, | float | Minimum longitude to extract from gridded data<sup>[1](#namelistfoot1)</sup> |
| lonmax | lonmax=360.0, | float | Maximum latitude to extract from gridded data<sup>[1](#namelistfoot1)</sup> |
| hrintvl | hrintvl=6.0, | float | Time interval (hours) between data points |
| trajinds | trajinds=0,1,2,3, | integer (4) | Indices in a TE track file corresponding to lon, lat, wind, pres |

<a name="namelistfoot1">1</a>: If using global data, -180. -> 360. allows for all lon orderings. For more regional definitions, this will be specific to the gridded dataset coordinate conventions.

There are two command line configurations.

```$> ncl ExTraTrack.ncl 'nlfile="nl.myconfig"' 'year_min_str="'${YYYY}'"' 'year_max_str="'${YYYY}'"'```

One, overrides the default namelist by specifying a user-defined namelist as `nlfile`. The second, `year_min_str` and `year_max_str` define the year boundaries for that particular instance of ExTraTrack. It can be used to parallelize the code (i.e., for a 20 year dataset you could spawn 20 single-core jobs doing 1984-1984, 1985-1985, and so on).

Simply, the code reads in the above TC trajectories and searches gridded data matching the period of analysis. During the times a TC trajectory exists, it calculates Hart's Vut, Vlt, and B.

Following the termination of a defined warm-core TC trajectory, the code continues following features based on sea level pressure minima at subsequent time slices until the cyclone dissipates or is unable to be tracked (leaves domain, end of gridded timeseries, etc.)

`reanalysis_et_cyclone_traj.ncl` will output new, storm-by-storm trajectory files in `./text_files/` with the prefix `tmp_TYPE_XXX`. Timing files are also output as `timing_TYPE_XXX` in `./timing_files/`

Each storm trajectory will look like...
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

### 5.) Concatenate, postprocess, and analyze statistics

**NOTE:** All NCL files beginning with `et_` call the default namelist but have an optional command line arg of `nlfile` as with ExTraTrack that allows for the user to specify a different namelist file.

In `./text_files/` the individual `tmp_XXX_NNN.txt` files need to be concatenated into a singlular traj file. This can be done quickly by invoking
```$> ncl et_concat_trajs.ncl```

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
