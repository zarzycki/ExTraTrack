;========================================================
load "../functions/mask_tc.ncl"
load "../functions/bearing.ncl"
load "../functions/et_track_plot.ncl"
load "../functions/et_track_checks.ncl"
load "../functions/et_phase_plot.ncl"
load "../functions/cps_calcs.ncl"
load "../functions/namelist_extract.ncl"
load "../functions/fcns-extract-storm.ncl"
load "../functions/getTrajectoriesGeneric.ncl"

;========================================================
; reads in TC trajectories and data
; follow storm through ET
; compute B, Vlt, and Vut
; plot track and plot phase space diagram
; write new ET trajectory file
;========================================================
begin
print(" ")

defaultnl="./namelist"
if (isdefined("nlfile")) then
  print("reading namelist from "+nlfile)
else if ((.not. isdefined("nlfile")) .and. (fileexists(defaultnl))) then
  print("using default namelist "+defaultnl)
  nlfile=defaultnl
else
  print("no namelist file defined, exiting...")
  print("feed in namelist file on the command line as in README")
  exit
end if
end if

type = namelist_extract( nlfile, "type", "string")
tfile = namelist_extract( nlfile, "tfile", "string")
filelist = namelist_extract( nlfile, "filelist", "string")
basin = namelist_extract( nlfile, "basin", "integer")
nllatmin = namelist_extract( nlfile, "latmin", "float")
nllatmax = namelist_extract( nlfile, "latmax", "float")
nllonmin = namelist_extract( nlfile, "lonmin", "float")
nllonmax = namelist_extract( nlfile, "lonmax", "float")
hrintvl = namelist_extract( nlfile, "hrintvl", "float")
trajinds = namelist_extract( nlfile, "trajinds", "integer")

; Load default settings (change at your own risk...)
nldefaults="./defaults/defaults.nl"
critRad = namelist_extract( nldefaults, "critRad", "float")
eps = namelist_extract( nldefaults, "eps", "float")
nlminlatetc = namelist_extract( nldefaults, "minlatetc", "float")
maskspacing = namelist_extract( nldefaults, "maskspacing", "float")
cpslevs = namelist_extract( nldefaults, "cpslevs", "float")
extractBox = namelist_extract( nldefaults, "extractBox", "float")
vt_calc_method = namelist_extract( nldefaults, "vt_calc_method", "string")
Zvar = namelist_extract( nldefaults, "Zvar", "string")
npastdays = namelist_extract( nldefaults, "npastdays", "float")

print(" ")

; CMZ setting to continue tracking PSL min past traj end
KEEPGOING=True
IGNORE_DISSIPATORS=True

; select year(s) and months(s) of traj files to look at.
; if not defined on the command line, you can specify by hardcoding here, useful for debugging
; else, will read in command line vars year_min_str, year_max_str and convert to ints
if (.not. isdefined("year_min_str") .or. .not. isdefined("year_max_str"))
  year_min = 1980
  year_max = 2022
  print("No CLI year_min_str / year_max_str,   using hardcoded: "+year_min+"   "+year_max)
else
  year_min=toint(year_min_str)
  year_max=toint(year_max_str)
  print("Specified year_min_str / year_max_str,   using: "+year_min+"   "+year_max)
end if
if (.not. isdefined("month_min_str") .or. .not. isdefined("month_max_str"))
  month_min = 1
  month_max = 12
  print("No CLI month_min_str / month_max_str, using hardcoded: "+month_min+"   "+month_max)
else
  month_min=toint(month_min_str)
  month_max=toint(month_max_str)
  print("Specified month_min_str / month_max_str, using: "+month_min+"   "+month_max)
end if

;========================================================
; Internal default calcs
;========================================================

; pre-calculate vars for degree/radian conversions
pi = 4.*atan(1.)
d2r = pi/180.
r2d = (1./d2r)

; debugging
plot_test=False
plot_phase=True
plot_track=False
B_print_nc=False

;========================================================
; Do some sanity checking and internal setting
;========================================================

if (vt_calc_method .eq. "regline") then
  interp_Z_to_CPS=False
else if (vt_calc_method .eq. "simple") then
  interp_Z_to_CPS=True
else
  print("invalid vt_calc_method: "+vt_calc_method)
  print("currently regline and simple are supported")
  exit
end if
end if
print("Using vt_calc_method: "+vt_calc_method)
print("So ExTraTrack has set interp_Z_to_CPS to "+interp_Z_to_CPS)

;========================================================
; Build lookup table from files listed in filelist
;========================================================

filel=asciiread(filelist,-1,"string")
timeArr=buildTimeLookup(filel)

;========================================================
; Get time units/calendar attributes from "file 0"
;========================================================

ftimetmp = addfile(filel(0),"r")

if (isatt(ftimetmp->time,"units")) then
  ftimeunits=ftimetmp->time@units
else
  print("no valid time units attribute on netcdf files, exiting...")
  exit
end if

if (isatt(ftimetmp->time,"calendar")) then
  ftimecal=ftimetmp->time@calendar
else
  print("no valid time calendar attribute on netcdf files, setting to standard...")
  ftimecal="standard"
end if

print("ftimeunits: "+ftimeunits)
print("ftimecal: "+ftimecal)
delete(ftimetmp)

;========================================================
; Read in TC trajectories
;========================================================

print("Loading TC trajectories...")

; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
nstorms_tot = stringtoint(systemfunc("grep 'start' "+tfile+" | wc -l"))
print("   Total number of storms = "+nstorms_tot)

ncols = -1                    ; assign number of columns
maxLines = -1
headerDelimStr = "start"
tdata = getTrajectoriesGeneric(tfile,ncols,maxLines,headerDelimStr,False,False)
tdata_dims = dimsizes(tdata)
tdata_nstorms = tdata_dims(0)
tdata_nvars = tdata_dims(1)
tdata_ntimes = tdata_dims(2)

; All arrays beginning t_??? are nstorms x ntimes array containing relevant
; TC trajectory information from TempestExtremes
; Use user-defined trajinds to set lon, lat, wind, pres
t_lon = tofloat(tdata(:,trajinds(0),:))
t_lat = tofloat(tdata(:,trajinds(1),:))
t_wind = tofloat(tdata(:,trajinds(2),:))
t_pres = tofloat(tdata(:,trajinds(3),:))
if (max(t_pres) .gt. 1100.) then
  print("WARNING: found pressures >1100. in track file, assuming Pa -> correcting to hPa")
  t_pres = t_pres/100.
end if

; Assume last 4 vars are YYYY MM DD HH
t_yyyy=  toint(tdata(:,tdata_nvars-4,:))
t_mm= toint(tdata(:,tdata_nvars-3,:))
t_dd= toint(tdata(:,tdata_nvars-2,:))
t_hh= toint(tdata(:,tdata_nvars-1,:))

; get start/end arrays that are length nstorms
s_dur   = dim_sum_n(where(.not.ismissing(t_lon),1,0),1)
s_year  = t_yyyy(:,0)
s_month = t_mm(:,0)
s_day   = t_dd(:,0)
s_hour  = t_hh(:,0)
s_id    = ispan(0,nstorms_tot-1,1)
e_year  = new(nstorms_tot,integer)
e_month = new(nstorms_tot,integer)
do nn = 0,nstorms_tot-1
  e_year(nn)  = t_yyyy(nn,s_dur(nn)-1)
  e_month(nn) = t_mm(nn,s_dur(nn)-1)
end do

delete([/tdata,nstorms_tot,headerDelimStr,maxLines,ncols/])

;========================================================
; Select storms within desired time period
;========================================================

; find indices where data is between given year/month min and max
ind_storm = ind(((s_year .eq. year_min .and. s_month .ge. month_min) .or.  \
                 (s_year .gt. year_min .and. s_year .le. year_max)) .and. \
                ((s_year .eq. year_max .and. s_month .le. month_max) .or.  \
                 (s_year .ge. year_min .and. s_year .lt. year_max)))

delete([/year_min,year_max,month_min,month_max,e_year,e_month/])

ntc = dimsizes(ind_storm)
print("   Number of TCs to be evaluated = "+ntc)
print("TC trajectories loaded.")
print(" ")

;========================================================
; Read in data
;========================================================

; loop over each storm
do jtc = 0,ntc-1
  itc = jtc + ind_storm(0)
  itcc = s_id(itc)
  print("Storm Number: "+itcc)
  print("  "+s_year(itc)+" "+s_month(itc)+" "+s_day(itc))

  nrows = s_dur(itc)

  ;========================================================
  ; Reset lat bounds to namelist and filter by hemisphere
  ;========================================================

  latmin=nllatmin
  latmax=nllatmax
  lonmin=nllonmin
  lonmax=nllonmax

  orilat=t_lat(itc,0)
  orilon=t_lon(itc,0)

  if (orilat .gt. 0) then
    print("NH storm, restricting lats")
    latmin=0.0
  else
    print("SH storm, restricting lats")
    latmax=0.0
  end if

  ;========================================================
  ; Test for correct basin
  ;========================================================

  if (.not. IGNORE_DISSIPATORS) then
    minlatetc=0.0
  else
    minlatetc=nlminlatetc
  end if

  ; select appropriate basin based on trajectory lat/lon
  ; analyze only north atlantic basin storms, skip all others

  if (basin .gt. 0) then
    print("  Checking storm basin...")

    tcMask = mask_tc(0,-999,-999,maskspacing,360)

    if (orilon .gt. 360.0-maskspacing) then
      orilon=360.0-maskspacing
      print("    original lon: "+t_lon(itc,0)+"   modified lon: "+orilon)
    end if

    if (tcMask({orilat},{orilon}) .ne. basin) then
      print("  user-specified basin: "+tcMask({orilat},{orilon})+" != this traj basin: "+basin+ \
        "... skipping storm")
      print(" ")
      delete([/nrows,tcMask,orilat,orilon/])
      continue
    else if (abs(t_lat(itc,s_dur(itc)-1)) .lt. minlatetc) then
      print("      Storm dissipates... skipping storm")
      print(" ")
      delete([/nrows,tcMask/])
      continue
    end if
    end if

    print("     "+tcMask({orilat},{orilon})+" = "+basin+"... matching basin, tracking storm")
    print(" ")

    delete([/tcMask/])

  end if

  STARTTIME=stringtointeger(systemfunc("date -u +%s"))

  ;========================================================
  ; Find start time
  ;========================================================

  ; convert longitude (0 to 360) to (-180 to 180)
  t_lon(itc,:) = where(t_lon(itc,:) .gt. 180, t_lon(itc,:) - 360, t_lon(itc,:))

  ; convert start time to match data calendar format

  OPT=0
  OPT@calendar := ftimecal
  s_time = cd_inv_calendar(s_year(itc),s_month(itc),s_day(itc),s_hour(itc),0,0, \
                           ftimeunits,OPT)
  s_date = sprinti("%0.4i",s_year)+sprinti("%0.2i",s_month)+sprinti("%0.2i",s_day)

  ; allow npastdays days past the end of TC phase
  if ( str_match_bool_ic(ftimeunits,"days") ) then
    e_time = s_time + (s_dur(itc)/4.) + npastdays
  else if ( str_match_bool_ic(ftimeunits,"hours") ) then
    e_time = s_time + (s_dur(itc)*6.) + (npastdays*24.)
  else
    print("confused about time units, check time units on input netcdf data, exiting...")
    exit
  end if
  end if
  e_time@units := s_time@units
  e_time@calendar := s_time@calendar

  ; convert time to date format for selecting files
  end_date = cd_calendar(e_time,0)
  e_mon_str = sprinti("%0.2i",floattoint(end_date(0,1)))
  e_year = end_date(0,0)

  print("TC start time: "+cd_calendar(s_time,-3))
  print("Extended end time: "+cd_calendar(e_time,-3))

  delete(end_date)

  ;========================================================
  ; Load input files
  ;========================================================

  if (s_time .lt. min(timeArr)) then
    print("Traj sttime: "+s_time+"  min gridded time: "+min(timeArr))
    print("start of traj lies outside available gridded data times, skipping...")
    delete([/nrows/])
    continue
  else
    if (e_time .gt. max(timeArr)) then
      print("Traj entime: "+e_time+"  max gridded time: "+max(timeArr))
      print("Setting end time to last time in array")
      print("WARNING: may not be full traj")
      e_time = max(timeArr)
    end if
    this_tmp1=getFileAndTimeFromLookup(timeArr,s_time)
    this_tmp2=getFileAndTimeFromLookup(timeArr,e_time)
    print("... files: "+this_tmp1[0]+" to "+this_tmp2[0])
    print("... file indices: "+ind(filel.eq.this_tmp1[0])+" to "+ind(filel.eq.this_tmp2[0]))
    files1=filel(ind(filel.eq.this_tmp1[0]):ind(filel.eq.this_tmp2[0]))
    delete([/this_tmp1,this_tmp2,e_year,e_mon_str/])
  end if

  ;========================================================
  ; Load files
  ;========================================================

  print("Loading files...")
	; read in data for Z levels
	input1 = addfiles(files1,"r")
	ListSetType(input1,"cat")
	delete(files1)
  print("... done loading files.")

  ;========================================================
  ; Figure out level name
  ;========================================================

  ; I believe this checks in descending priority, but I can't imagine
  ; a way for a Z to have multiple vertical dimensions?
  if isfilevarcoord( input1[0], Zvar, "plev") then
    levname="plev"
  else if isfilevarcoord( input1[0], Zvar, "lev_p") then
    levname="lev_p"
  else if isfilevarcoord( input1[0], Zvar, "lev") then
    levname="lev"
  else if isfilevarcoord( input1[0], Zvar, "level") then
    levname="level"
  else
    print("Unsupported level dimension on "+Zvar+" file")
    print("please name it lev or level")
    exit
  end if
  end if
  end if
  end if
  print("levname: "+levname)

  ;========================================================
  ; Check Z levs
  ;========================================================

  ; Check to make sure levels are oriented from top to bottom (i.e., low P to high P)
  flipLev=False
  lev_in = input1[0]->$levname$
  lev = tofloat(lev_in)
 	nlev = dimsizes(lev)
  if (lev(0) .gt. lev(nlev-1)) then
    print("WARNING: levs ordered from bottom (1000mb) to top (300mb)")
    print("WARNING: attempting to reverse, but you should check this!")
    flipLev=True
  end if
	delete([/lev,lev_in/])

  ;========================================================
  ; Load Z
  ;========================================================

	print("Loading geopotential heights")
  if (flipLev) then
    Z3 = input1[:]->$Zvar$({s_time:e_time},::-1,{latmin:latmax},:)
	  ;Z3 = input1[:]->Z({s_time:e_time},{latmin:latmax},:,::-1)
  else
    Z3 = input1[:]->$Zvar$({s_time:e_time},:,{latmin:latmax},:)
	  ;Z3 = input1[:]->Z({s_time:e_time},{latmin:latmax},:,:)
  end if

  ; make sure we are ordered correctly
  Z3 := Z3(time|:,$levname$|:,lat|:,lon|:)

	if (max(Z3) .gt. 20000.) then
	  print("WARNING: Z seems high, do you have geopotential instead of geo height?")
	  print("WARNING: dividing by g, but you should check this!")
	  Z3 = Z3/9.81
	end if
	nlev = dimsizes(Z3&$levname$)

	if (interp_Z_to_CPS) then
	  print("interpolating Z to CPS surfaces")
	  ; convert cpslevs to same type as Z3 levs
	  cpslevs := totype(cpslevs,typeof(Z3&$levname$))
	  Z3 := int2p_n_Wrap(Z3&$levname$,Z3,cpslevs,1,1)
	  Z3!1=levname
	end if

  ; pressure levels
  pnew = Z3&$levname$
  logp = log(pnew*100.)   ; convert p from hPa to Pa and take log

  ; assign indices to be used later
	indtop = closest_val(cpslevs(0),pnew)      ; top index for VTU
	indmid = closest_val(cpslevs(1),pnew)      ; bottom index for VTU, top index for B, VTL
	indbot = closest_val(cpslevs(2),pnew)      ; bottom index for B, VTL
  print("CPS indices: indtop: "+indtop+"  indmid: "+indmid+"  indbot: "+indbot)

	; compute weights from latitudes
	wgt = cos((Z3&lat)*d2r)

  ; flip longitude for geopotential height
	Z = Z3(:,:,:,{lonmin-eps:lonmax})
  delete([/Z3/])
  print("... done loading geopotential height")
  ;========================================================
  ; Load data and calculate UMAG
  ;========================================================

  print("Loading U+V...")
	; load wind speeds at lowermost model level
	UBOT_in = input1[:]->UBOT({s_time:e_time},{latmin:latmax},:)
	VBOT_in = input1[:]->UBOT({s_time:e_time},{latmin:latmax},:)

  ; flip longitude
	UBOT = UBOT_in(:,:,{lonmin-eps:lonmax})
	VBOT = VBOT_in(:,:,{lonmin-eps:lonmax})
	delete([/UBOT_in,VBOT_in/])

  UBOTMAG = ((VBOT^2 + UBOT^2)^(0.5))

  copy_VarMeta(UBOT,UBOTMAG)
	delete([/UBOT,VBOT/])
  print("... done loading U+V.")

  ;========================================================
  ; Load data and preallocate some arrays
  ;========================================================

  print("Loading rest of data...")
  ; time and lat/lon arrays
  time = input1[:]->time({s_time:e_time})
  ; because lat/lon may be double and we want float, but NCL doesn't persist meta when
  ; type converting, we load temp arrays, convert to float, and copy meta, release tmp arrays
  lat_in = input1[0]->lat({latmin:latmax})
  lon_in = input1[0]->lon({lonmin-eps:lonmax})
  lat = tofloat(lat_in)
  lon = tofloat(lon_in)

  ; attempt to check whether data is -180-180 or 0-360.
  ; lonpos = True  --> 0-360
  ; lonpos = False --> -180-180
  ; NOTE: for regional data, we may not know, set to False to keep lat/lon as in traj file
  if (max(lon) .gt. 180.0) then
    lonpos=True
  else if (min(lon) .lt. 0.0) then
    lonpos=False
  else
    lonpos=False
  end if
  end if
  copy_VarMeta(lat_in,lat)
  copy_VarMeta(lon_in,lon)
	delete([/lat_in,lon_in/])

	; flip longitude for sea level pressure
	PSL_all = input1[:]->PSL({s_time:e_time},{latmin:latmax},:)
	PSL_in = PSL_all(:,:,{lonmin-eps:lonmax})

  date = cd_calendar(input1[:]->time({s_time:e_time}),-2)
	delete([/PSL_all,input1,s_time,e_time/])

  ; convert Pa to hPa
  PSL = PSL_in/100.0
  copy_VarMeta(PSL_in,PSL)
  PSL@units := "hPa"
  delete(PSL_in)

  ; assign fill value for missing if doesn't exist
  if(.not.isatt(PSL,"_FillValue")) then
    print("assigning fillvalue to PSL")
    PSL@_FillValue := 9.969209968386869e+36
  end if
  if(.not.isatt(Z,"_FillValue")) then
    print("assigning fillvalue to Z")
    Z@_FillValue := 9.969209968386869e+36
  end if
  if(.not.isatt(UBOTMAG,"_FillValue")) then
    print("assigning fillvalue to UBOTMAG")
    UBOTMAG@_FillValue := 9.969209968386869e+36
  end if
  print("... done loading rest of data")

  ; get Z thickness between 600 and 900 for B calculation
  Z_thick = Z(:,indmid,:,:) - Z(:,indbot,:,:)
  Z_thick!0 = "time"
  Z_thick!1 = "lat"
  Z_thick!2 = "lon"
  Z_thick&time = time
  Z_thick&lat = lat
  Z_thick&lon = lon

  ;========================================================
  ; Preallocate arrays
  ;========================================================

  ; record dimension sizes
  ntime = dimsizes(time)  ; time
  nlat = dimsizes(lat)    ; latitude
  nlon = dimsizes(lon)    ; longitude

  ; preallocate arrays for minimum pressure value and locations
  ; note that ntime corresponds to the netCDF data loaded *not* the TE trajectory files
  minP    = new(ntime,float)   ; minimum pressure value (Pa)
  minPlat = new(ntime,float)   ; minimum pressure latitude
  minPlon = new(ntime,float)   ; minimum pressure longitude
  V       = new(ntime,float)   ; maximum wind speed at 10 m
  dist    = new(ntime,float)   ; distance traveled between each time step (m)
  ang     = new(ntime,float)   ; angle of storm travel
  d_ang   = new(ntime,float)   ; change in angle of storm travel
  B       = new(ntime,float)   ; B parameter (Hart 2003)
  Vlt     = new(ntime,float)   ; lower troposphere thermal wind (Hart 2003)
  Vut     = new(ntime,float)   ; upper troposphere thermal wind (Hart 2003)
  ang_all = new((/nlat,nlon/),float) ; angle of all values from storm center

  ;========================================================
  ; Compute values during TC tracker
  ;========================================================

  print("Tracking storm...")

  ; set values to trajectory
  minP(0:nrows-1) = t_pres(itc,0:nrows-1)      ; minimum pressure values
  minPlat(0:nrows-1) = t_lat(itc,0:nrows-1)    ; latitude of minimum pressure
  minPlon(0:nrows-1) = t_lon(itc,0:nrows-1)    ; longitude of minimum pressure
  V(0:nrows-1) = t_wind(itc,0:nrows-1)         ; max wind speed
  
  ; check to see if we are missing any data
  OPT:=0
  OPT@calendar:=ftimecal
  ; dummy zero arrays for min + sec in cd_inv_calendar of same size as yyyy, mm, dd, hh
  t_zero = t_yyyy(itc,0:nrows-1)
  t_zero = 0
  ; use hours here because hrintvl is reported in hours.
  ; the exact "units" don't matter since we just want some float string that allows us to fill in
  ; gaps with equally spaced hours
  te_time = cd_inv_calendar(t_yyyy(itc,0:nrows-1),t_mm(itc,0:nrows-1),t_dd(itc,0:nrows-1),t_hh(itc,0:nrows-1),t_zero,t_zero,"hours since 0001-01-01 00:00",OPT)
  ; use max/min te_time to figure out how many data points *should* lie between these
  nints=toint( (max(te_time)-min(te_time)) / hrintvl ) + 1
  
  if (nints .gt. nrows) then
    print("INTERP: Found disagreement between TE trajectory and loaded NetCDF times")
    print("INTERP: This is typically due to TE having an allowable maxgap")
    print("INTERP: TE times: "+dimsizes(te_time)+"    interpolated times: "+nints)
    print("INTERP: Will attempt to linearly interpolate between "+min(te_time)+" and "+max(te_time))
    int_time = fspan(min(te_time),max(te_time),nints)
    nirows = dimsizes(int_time)
    minP(0:nirows-1)    = linint1(te_time,minP(0:nrows-1),False,int_time, 0)
    minPlat(0:nirows-1) = linint1(te_time,minPlat(0:nrows-1),False,int_time, 0)
    minPlon(0:nirows-1) = linint1(te_time,minPlon(0:nrows-1),False,int_time, 0)
    V(0:nirows-1)       = linint1(te_time,V(0:nrows-1),False,int_time, 0)
    ; recast nirows as nrows now that we've updated all arrays
    nrows = nirows
    delete([/int_time,nirows/])
  end if

  ; clean up
  delete([/t_zero,te_time,nints,OPT/])

  ; if we have 0->360 gridded data, let's check to see if we have -180-180 traj lon
  ; if so, change traj data to match gridded data coords
  if (lonpos .eq. True) then
    minPlon=where(minPlon.lt.0,360+minPlon,minPlon)
  end if

  do i = 1,nrows-1
    ;========================================================
    ; compute distance traveled, speed of storm, and change in direction of motion
    ; for times already in TC trajectory
    ;========================================================
    dist(i) = gc_latlon(minPlat(i),minPlon(i),minPlat(i-1),minPlon(i-1),2,4)  ; km
    ang(i) = bearing(minPlat(i-1),minPlon(i-1),minPlat(i),minPlon(i))
    if (i .eq. 1) then
      ; if first datapoint, we can't have a change in angle, so set to zero
      d_ang(i) = 0.
    else
      d_ang(i) = r2d*atan2(sin((ang(i-1)-ang(i))*d2r),cos((ang(i-1)-ang(i))*d2r))
    end if
    ;========================================================
    ; Set options to pass into CPS calcs
    ;========================================================
    opt=0
    opt@date:=cd_calendar(time(i),3)
    opt@id:=itc
    opt@v_interp=vt_calc_method  ; regline or simple
    if (B_print_nc) then
      opt@debug=True
    else
      opt@debug=False
    end if

    ;========================================================
    ; Find thickness of the right and left hand sides of storm
    ;========================================================
    B(i) = calc_B(Z_thick(i,:,:),extractBox,critRad,dist(i),ang(i),minPlat(i),minPlon(i),opt)

    ;========================================================
    ; Compute thermal wind values (Vlt and Vut)
    ;========================================================
    rtn_tmp = calc_VltVut(Z(i,:,:,:),lat,lon,extractBox,critRad,dist(i), \
      ang(i),minPlat(i),minPlon(i),(/indtop,indmid,indbot/),logp,opt)
    Vlt(i) = rtn_tmp(0)
    Vut(i) = rtn_tmp(1)
    delete(rtn_tmp)

    ;========================================================
    ; Print various fields to screen
    ;========================================================
    if (i .eq. 1) then
      print("    Press.   Lat      Lon     dist    V      ang" + \
            "     d_ang      B     Vlt     Vut     TYPE")
    end if
      print_ET_toScreen(minP(i),minPlat(i),minPlon(i),dist(i),V(i),ang(i),d_ang(i), \
        B(i),Vlt(i),Vut(i),"TC")

  end do   ; end loop of already TC tracked stuff

  ;========================================================
  ; call function to create et plot
  ;========================================================

  if (plot_track) then
    do ii = 0, nrows-1
      et_track_plot(type,itc,ii,lat,lon,date(ii),minP(ii),minPlat(ii),minPlon(ii), \
                    PSL(ii,:,:),PSL(ii,:,:))
    end do
    delete(ii)
  end if

  if (KEEPGOING) then

  ;========================================================
  ; Find great circle distances and pressure values within radius
  ;========================================================

  ; loop over all times after TC tracker ends
  do i = nrows,ntime-1
    if (minPlat(i-1) .ge. (latmax-10.) .or. minPlon(i-1) .ge. (lonmax-10.) .or. \
           minPlat(i-1) .le. (latmin+10.) .or. minPlon(i-1) .le. (lonmin+10.)) then
      print("Done tracking storm since it exited our available domain...")
      break
    else  ; still tracking storm
      ; find distance from storm center to each grid point using great circle length
      PSL_tmp = PSL(i,:,:)
      UBOTMAG_tmp = UBOTMAG(i,:,:)

      ; get radius for storm size given latitude
      ; Mask PS and UBOT arrays based on radius as a function of latitude.
      rad_tmp=get_sizeRad_from_lat(minPlat(i-1))
      lat2d=conform_dims((/nlat,nlon/),lat,0)
      lon2d=conform_dims((/nlat,nlon/),lon,1)
      d = gc_latlon(minPlat(i-1),minPlon(i-1),lat2d,lon2d,2,4)  ; km
      PSL_tmp=where(d.gt.rad_tmp,PSL_tmp@_FillValue,PSL_tmp)
      UBOTMAG_tmp=where(d.gt.rad_tmp,UBOTMAG_tmp@_FillValue,UBOTMAG_tmp)

      ; find max wind speed
      V_tmp = max(UBOTMAG_tmp)

      ; clean up
      delete([/d,lat2d,lon2d,rad_tmp,UBOTMAG_tmp/])

      ;========================================================
      ; Find all local pressure minimums
      ;========================================================

      ; find location of pressure minimums
      locmin = local_min(PSL_tmp,False,0.0)

      ; if pressure minimums exist, choose the best one
      if (locmin(0) .gt. 0) then
        locmin_lat = lat(locmin@yi)    ; latitudes of local mins
        locmin_lon = lon(locmin@xi)    ; longitudes of local mins
        nmins = dimsizes(locmin_lat)   ; number of local mins

        ; plot all local mins for debugging
        plot_test = False
        if (plot_test) then
          plot_track = False   ; don't overwrite with later plot of single point
          et_track_plot(type,itc,i,lat,lon,date(i),minP(i-1),locmin_lat,locmin_lon, \
                PSL_tmp,PSL(i,:,:))
        end if

        ;========================================================
        ; Test if each local min is a valid storm
        ;========================================================

        ; preallocate arrays
        locmin_test = new(nmins,logical)  ; test if local min is valid
        minP_tmp = new(nmins,float)       ; minimum pressure value
        dist_tmp = new(nmins,float)       ; distance traveled (km)
        ang_tmp = new(nmins,float)        ; storm bearing (degrees)
        d_ang_tmp = new(nmins,float)      ; change in bearing (degrees)



        ; loop over all local minimums
        do jmin = 0,nmins-1

           minP_tmp(jmin) = PSL(i,{locmin_lat(jmin)},{locmin_lon(jmin)})
           dist_tmp(jmin) = gc_latlon(minPlat(i-1),minPlon(i-1), \
                     locmin_lat(jmin),locmin_lon(jmin),2,4) ; km
           ang_tmp(jmin) = bearing(minPlat(i-1),minPlon(i-1), \
                       locmin_lat(jmin),locmin_lon(jmin))
           d_ang_tmp(jmin) = r2d*atan2(sin((ang(i-1)-ang_tmp(jmin))*d2r), \
                         cos((ang(i-1)-ang_tmp(jmin))*d2r))
           locmin_test(jmin) = et_track_checks(latmax,lonmax, \
                               minP_tmp(jmin),locmin_lat(jmin),locmin_lon(jmin), \
                     dist_tmp(jmin),d_ang_tmp(jmin))
        end do

        delete(jmin)

        locmin_true = ind(locmin_test)   ; gather indices where storm is valid


        ;========================================================
        ; Choose best storm (valid + lowest pressure)
        ;========================================================

        if (.not. all(ismissing(locmin_true))) then
          minPind = ind(minP_tmp .eq. min(minP_tmp(locmin_true)))  ; choose lowest pressure
          ; loop to check if two valid points within the search area
          ; have the exact same minimum pressure (yes, this actually happens)
          ; choose minimum pressure with travel distance/speed most similar to previous time step
          if (dimsizes(minPind) .gt. 1) then
            dist_diff = abs(dist(i-1) - dist_tmp)
            minPindTmp=-1
            minDiffTmp=99999.9
            do ii = 0,dimsizes(minPind)-1
              print(dist_diff(minPind(ii))+" "+minDiffTmp)
              if (dist_diff(minPind(ii)) .lt. minDiffTmp) then
                minDiffTmp=dist_diff(minPind(ii))
                print("updating minPindTmp from "+minPindTmp+" "+minPind(ii))
                minPindTmp=minPind(ii)
              end if
            end do
            delete([/minDiffTmp,dist_diff,ii,minPind/])
            minPind=minPindTmp
            delete(minPindTmp)
          end if

          minPlat(i) = (/locmin_lat(minPind)/)   ; latitude of minimum pressure
          minPlon(i) = (/locmin_lon(minPind)/)   ; longitude of minimum pressure
          minP(i) = minP_tmp(minPind)            ; minimum pressure value
          dist(i) = dist_tmp(minPind)     ; distance traveled by storm (km)
          V(i) = V_tmp                    ; max wind speed
          ang(i) = ang_tmp(minPind)       ; storm bearing
          d_ang(i) = d_ang_tmp(minPind)   ; change in storm direction

          delete(minPind)
        else ; stop tracking if no local min satisfies all tests
          delete([/PSL_tmp,locmin,locmin_lat,locmin_lon,nmins/])
          delete([/locmin_test,minP_tmp,dist_tmp,V_tmp,ang_tmp,d_ang_tmp,locmin_true/])
          break
        end if

        delete([/locmin,locmin_lat,locmin_lon,nmins/])
        delete([/locmin_test,minP_tmp,dist_tmp,V_tmp,ang_tmp,d_ang_tmp,locmin_true/])

      else  ; stop tracking if no local min exists
        delete([/PSL_tmp,locmin/])
        break
      end if

      ;========================================================
      ; Set options to pass into CPS calcs
      ;========================================================
      opt=0
      opt@date:=cd_calendar(time(i),3)
      opt@id:=itc
      if (B_print_nc) then
        opt@debug:=True
      else
        opt@debug:=False
      end if

      ;========================================================
      ; Find thickness of the right and left hand sides of storm
      ;========================================================
      B(i) = calc_B(Z_thick(i,:,:),extractBox,critRad,dist(i),ang(i),minPlat(i),minPlon(i),opt)

      ;========================================================
      ; Compute thermal wind values (Vlt and Vut)
      ;========================================================
      rtn_tmp = calc_VltVut(Z(i,:,:,:),lat,lon,extractBox,critRad,dist(i), \
        ang(i),minPlat(i),minPlon(i),(/indtop,indmid,indbot/),logp,opt)
      Vlt(i) = rtn_tmp(0)
      Vut(i) = rtn_tmp(1)
      delete(rtn_tmp)

      ;========================================================
      ; Print various fields to screen
      ;========================================================
      print_ET_toScreen(minP(i),minPlat(i),minPlon(i),dist(i), \
        V(i),ang(i),d_ang(i),B(i),Vlt(i),Vut(i),"ETC")

      ;========================================================
      ; Plot storm
      ;========================================================

      ; call function to create et plot
      if (plot_track) then
        system("mkdir -p ./gifs")
        et_track_plot(type,itc,i,lat,lon,date(i), \
                      minP(i),minPlat(i),minPlon(i),PSL_tmp,PSL(i,:,:))
      end if

      delete(PSL_tmp)

    end if ; check to see if storm is still in domain
  end do   ;loop over all times after TC tracker ends

  end if  ; end KEEPGOING

  delete([/PSL,ang_all,Z,Z_thick,UBOTMAG,wgt,lat,lon,nlat,nlon,ntime,nrows/])

  print("Finished tracking storm.")
  print(" ")

  ; if plotting track, combine pdfs into gif
  if (plot_track .or. plot_test) then
    system("mkdir -p ./gifs")
    system("rm -f ./gifs/track_plot_"+type+"_"+sprinti("%0.4i",itcc)+".gif")
    system("convert -delay 50 -loop 0 -density 300 ./gifs/track_plot_"+type+"_"+ \
           sprinti("%0.4i",itcc)+"_*.pdf ./gifs/track_plot_"+type+"_"+ \
           sprinti("%0.4i",itcc)+".gif")
    system("rm -f ./gifs/track_plot_"+type+"_"+sprinti("%0.4i",itcc)+"_*.pdf")
  end if

  ;========================================================
  ; Plot phase space diagram
  ;========================================================

  if (plot_phase) then
    phaseplotdir="./phase_plots/"
    if ( toint(systemfunc("test -d "+phaseplotdir+"; echo $?")) .ne. 0 ) then
      system("mkdir -p "+phaseplotdir+"/")
    end if
    system("rm -f "+phaseplotdir+"/phase_plot_"+type+"_"+itc+".pdf") ; remove existing pdf
    et_phase_plot(type,itc,date(0:i-1),B,Vlt,Vut)     ; plot phase space diagram
  end if

  delete([/date/])

  ;========================================================
  ; Write output to text file
  ;========================================================

  print("Writing to trajectory ascii file...")

  textfiledir="./text_files/"
  timingfiledir="./timing_files/"

  ; make text and timing file folders if they do not exist.
  if ( toint(systemfunc("test -d "+textfiledir+"; echo $?")) .ne. 0 ) then
    system("mkdir -p "+textfiledir+"/")
  end if
  if ( toint(systemfunc("test -d "+timingfiledir+"; echo $?")) .ne. 0 ) then
    system("mkdir -p "+timingfiledir+"/")
  end if

  ; remove text file if it already exists
  system("rm -f "+textfiledir+"/tmp_"+type+"_"+sprinti("%0.4i",itcc)+".txt")

  ; write header line for storm
  write_table(textfiledir+"/tmp_"+type+"_"+sprinti("%0.4i",itcc)+".txt","w", \
              [/"start",i,s_year(itc),s_month(itc),s_day(itc),s_hour(itc),itcc/], \
              "%s   %03i  %04i      %02i    %02i    %02i    %04i")

  ; reformat missing values to avoid terrible formatting in ascii file
  dist(0) = -999
  ang(0)  = -999
  B(0)    = -999
  Vlt(0)  = -999
  Vut(0)  = -999

  st_date = cd_calendar(time(0:i-1),0)  ; get date at each time step

  ; write storm data to ascii file (append to header)
  write_table(textfiledir+"/tmp_"+type+"_"+sprinti("%0.4i",itcc)+".txt","a", \
    [/minPlon(0:i-1),minPlat(0:i-1),minP(0:i-1), \
    V(0:i-1),dist(0:i-1),ang(0:i-1), \
    B(0:i-1),Vlt(0:i-1),Vut(0:i-1), \
    st_date(:,0),st_date(:,1),st_date(:,2),st_date(:,3)/], \
    "   %6.2f   %5.2f   %7.2f" + \
    "   %4.1f   %7.2f   %7.2f" + \
    "   %7.2f   %7.2f   %7.2f" + \
    "   %4.0f  %2.0f  %2.0f  %2.0f")

  ENDTIME=stringtointeger(systemfunc("date -u +%s"))
  TOTTIME=ENDTIME-STARTTIME
  asciiwrite (timingfiledir+"/timing_"+type+"_"+sprinti("%0.4i",itcc)+".TIME",TOTTIME)

  delete([/i,itc,itcc,time,st_date,minPlon,minPlat,minP,V,dist,ang,d_ang,B,Vlt,Vut/])

  print("Done with storm.")
  print(" ")

end do ; end storm loop

end
