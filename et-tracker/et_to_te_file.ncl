;========================================================
load "../functions/mask_tc.ncl"
load "../functions/namelist_extract.ncl"
load "../functions/getTrajectoriesGeneric.ncl"
load "../functions/getStormIndices.ncl"
;========================================================
; USAGE:
;ncl et_to_te_file.ncl 'nlfile="user-nl/nl.hyp.ref.dtime900.002"'
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

etfile = namelist_extract( nlfile, "etfileavg", "string")

te_outfile = str_sub_str(etfile,"avg","te")
gridfile="/glade/u/home/zarzycki/work/unigridFiles/ne0np4natlanticref.ne30x4/topo/topo_ne0np4natlanticref.ne30x4_smooth.nc"

print("Loading TC trajectories")

; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
nstorms_tot = stringtoint(systemfunc("grep 'start' "+etfile+" | wc -l"))
print("   Total number of storms = "+nstorms_tot)

ncols = -1                    ; assign number of columns
maxLines = -1
headerDelimStr = "start"
tdata = getTrajectoriesGeneric(etfile,ncols,maxLines,headerDelimStr,False,False)
tdata_dims = dimsizes(tdata)
tdata_nstorms = tdata_dims(0)
tdata_nvars = tdata_dims(1)
tdata_ntimes = tdata_dims(2)

; All arrays beginning t_??? are nstorms x ntimes array containing relevant
; TC trajectory information from TempestExtremes
; Use user-defined trajinds to set lon, lat, wind, pres
trajinds=(/0,1,2,4/)
t_lon = tofloat(tdata(:,trajinds(0),:))
t_lat = tofloat(tdata(:,trajinds(1),:))

; Assume last 4 vars are YYYY MM DD HH
t_yyyy=  toint(tdata(:,tdata_nvars-4,:))
t_mm= toint(tdata(:,tdata_nvars-3,:))
t_dd= toint(tdata(:,tdata_nvars-2,:))
t_hh= toint(tdata(:,tdata_nvars-1,:))

wind  = tofloat(tdata(:,3,:))
psl   = tofloat(tdata(:,2,:))
dist  = tofloat(tdata(:,4,:))
ang   = tofloat(tdata(:,5,:))
B     = tofloat(tdata(:,6,:))
Vlt   = tofloat(tdata(:,7,:))
Vut   = tofloat(tdata(:,8,:))

; set fill values
fillvalue = -999.00
dist@_FillValue = fillvalue
ang@_FillValue = fillvalue
B@_FillValue = fillvalue
Vlt@_FillValue = fillvalue
Vut@_FillValue = fillvalue

if (max(psl) .gt. 1100.) then
  print("WARNING: found pressures >1100. in track file, assuming Pa -> correcting to hPa")
  psl = psl/100.
end if

; get start/end arrays that are length nstorms
s_dur   = dim_sum_n(where(.not.ismissing(t_lon),1,0),1)
s_year  = t_yyyy(:,0)
s_month = t_mm(:,0)
s_day   = t_dd(:,0)
s_hour  = t_hh(:,0)
s_id    = ispan(0,nstorms_tot-1,1)
e_year  = new(nstorms_tot,integer)
e_month = new(nstorms_tot,integer)
do ii = 0,nstorms_tot-1
  e_year(ii)  = t_yyyy(ii,s_dur(ii)-1)
  e_month(ii) = t_mm(ii,s_dur(ii)-1)
end do

; get storm indices from etfile header rows
s_num = getStormIndices(etfile)


; Write a TE-compatible file

quote = inttochar(34)
delim="	"
system("rm -rf "+te_outfile)

if (gridfile .ne. "") then
  gridf=addfile(gridfile,"r")
  gridlat=gridf->lat
  gridlon=gridf->lon
else
  print("no gridfile")
  exit
end if

print_to_screen=True
header="start"

do ii = 0,tdata_nstorms-1

  HEADSTR=header+delim+s_dur(ii)+delim+s_year(ii)+delim+s_month(ii)+delim+s_day(ii)+delim+s_hour(ii)

  if print_to_screen then
    print(HEADSTR+"")
  end if
  system("echo "+quote+      HEADSTR      +quote+" >> "+te_outfile)

  do jj = 0,s_dur(ii)-1
    ; if gridfile not empty string, find nearest coordinate
    if (gridfile .ne. "") then
      gcdist = gc_latlon(t_lat(ii,jj),t_lon(ii,jj),gridlat,gridlon,2,4)
      thisLat=minind(gcdist)
      thisLon=1
    else
      thisLat=-999
      thisLon=-999
    end if

    STORMSTR=delim+thisLat+delim+sprintf("%6.2f", t_lon(ii,jj))+delim+sprintf("%6.2f", t_lat(ii,jj))+delim+t_yyyy(ii,jj)+delim+t_mm(ii,jj)+delim+t_dd(ii,jj)+delim+t_hh(ii,jj)

    if print_to_screen then
       print(STORMSTR+"")
    end if

    system("echo "+quote+      STORMSTR       +quote+" >> "+te_outfile)
  end do
end do




end
