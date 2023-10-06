;========================================================
load "../functions/namelist_extract.ncl"
load "../functions/getTrajectoriesGeneric.ncl"
load "../functions/getStormIndices.ncl"
;========================================================
; reads in ET trajectories from "raw" tracker
; performs a running mean on B, Vut, and Vlt fields
; to smooth phase space
;
; orig code: D. Thatcher 2014/08/11
; v0.1 release/cleanup: C. Zarzycki 2017/08/08
;
;========================================================
begin
print("et_avg_text.ncl")
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

;========================================================
; User options
;========================================================

; unsmoothed ET traj file
etfile = namelist_extract( nlfile, "etfileori", "string")

; type string (used for defining output), generally a model shortname
typestr = namelist_extract( nlfile, "type", "string")

; frequency of output data (hrs)
datafreq = namelist_extract( nlfile, "hrintvl", "float")

; average file
catfile = namelist_extract( nlfile, "etfileavg", "string")

;========================================================
; Default namelist settings
;========================================================

nldefaults="./defaults/defaults.nl"

; running average time period (hrs)
hourmean = namelist_extract( nldefaults, "hourmean", "float")

;**********************************************************************************

;========================================================
; Read in ET storm trajectories
;========================================================

;========================================================
; Read in TC-only trajectories from traj file
;========================================================

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

year_min=min(s_year)
year_max=max(s_year)

delete([/tdata,headerDelimStr,maxLines,ncols,ii/])

;========================================================
; Loop over each storm
;========================================================

do jj = 0,nstorms_tot-1

  ; itc is the stormID (i.e., corresponds to integer index from original TE file)
  ; jtc is the index of the matching stormID in the etfile
  ; NOTE: unless we are filtering here, jtc should = i
  itc = s_num(jj)
  jtc = ind(s_num .eq. itc)

  print("Storm Number: "+itc+" with jtc: "+jtc)

  ; get duration of storm so arrays are only valid values and not missing
  dur = s_dur(jtc)

  ; Extract this storm's info from the master storm arrays
  lon_sub  = t_lon(jtc,0:dur-1)
  lat_sub  = t_lat(jtc,0:dur-1)
  pres_sub = psl(jtc,0:dur-1)
  V_sub    = wind(jtc,0:dur-1)
  dist_sub = dist(jtc,0:dur-1)
  ang_sub  = ang(jtc,0:dur-1)
  B_sub    = B(jtc,0:dur-1)
  Vlt_sub  = Vlt(jtc,0:dur-1)
  Vut_sub  = Vut(jtc,0:dur-1)
  ; NOTE: we are going to convert this back to float because of Diana's
  ; original formatting below and the desire to not change that.
  year_sub = tofloat(t_yyyy(jtc,0:dur-1))
  mon_sub  = tofloat(t_mm(jtc,0:dur-1))
  day_sub  = tofloat(t_dd(jtc,0:dur-1))
  hour_sub = tofloat(t_hh(jtc,0:dur-1))

  ; set fill values
  fillvalue = -999.00
  B_sub@_FillValue = fillvalue
  Vlt_sub@_FillValue = fillvalue
  Vut_sub@_FillValue = fillvalue

  ; allow for missing data <= -990.0
  B_sub = where(B_sub .le. -990.00,B_sub@_FillValue,B_sub)
  Vlt_sub = where(Vlt_sub .le. -990.00,Vlt_sub@_FillValue,Vlt_sub)
  Vut_sub = where(Vut_sub .le. -990.00,Vut_sub@_FillValue,Vut_sub)

  ;========================================================
  ; running avg
  ;========================================================

  ; perform running avg over hourmean hours
  nave = toint(hourmean/datafreq) + 1
  opt = 1    ; reflective (symmetric) conditions

  ; First, set all variables to un-averaged arrays
  B_avg = B_sub
  Vlt_avg = Vlt_sub
  Vut_avg = Vut_sub

  ; Now average excluding the 0th index, which has bad data
  if (dur .gt. nave) then
    print("... averaging")
  else
    print("... reducing average from "+nave+" to "+(dur - 1))
    nave = dur - 1
  end if

  if (dur .ge. 2)
    B_avg(1::) = runave_Wrap(B_sub(1::),nave,opt)
    Vlt_avg(1::) = runave_Wrap(Vlt_sub(1::),nave,opt)
    Vut_avg(1::) = runave_Wrap(Vut_sub(1::),nave,opt)
  end if

  delete([/nave,opt,B_sub,Vlt_sub,Vut_sub/])

  ;========================================================
  ; Write output to text file
  ;========================================================

  print("Writing to trajectory ascii file...")

  outfile = "tmp_avg_"+typestr+"_"

  textfiledir="./text_files/"
  ; make text and timing file folders if they do not exist.
  if ( toint(systemfunc("test -d "+textfiledir+"; echo $?")) .ne. 0 ) then
    system("mkdir "+textfiledir+"/")
  end if

  ; remove text file if it already exists
  system("rm -f "+textfiledir+"/"+outfile+sprinti("%0.4i",itc)+".txt")

  ; write header line for storm
  write_table(textfiledir+"/"+outfile+sprinti("%0.4i",itc)+".txt","w", \
              [/"start",dur,year_sub(0),mon_sub(0),day_sub(0),hour_sub(0),itc/], \
              "%s   %03i  %4.0f      %02.0f    %02.0f    %02.0f    %04i")

  ; write storm data to ascii file (append to header)
  write_table(textfiledir+"/"+outfile+sprinti("%0.4i",itc)+".txt","a", \
    [/lon_sub,lat_sub,pres_sub, \
    V_sub,dist_sub,ang_sub, \
    B_avg,Vlt_avg,Vut_avg, \
    year_sub,mon_sub,day_sub,hour_sub/], \
    "   %6.2f   %5.2f   %7.2f" + \
    "   %4.1f   %7.2f   %7.2f" + \
    "   %7.2f   %7.2f   %7.2f" + \
    "   %4.0f  %02.0f  %02.0f  %02.0f")

  delete([/dur,itc,lon_sub,lat_sub,pres_sub,V_sub,dist_sub,ang_sub/])
  delete([/B_avg,Vlt_avg,Vut_avg,year_sub,mon_sub,day_sub,hour_sub/])

  print("Done with storm.")
  print(" ")

end do  ; end storm loop

system("cat "+textfiledir+"/"+outfile+"* > "+catfile)
system("rm -f "+textfiledir+"/"+outfile+"*")

end
