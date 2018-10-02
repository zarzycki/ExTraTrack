;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "../functions/namelist_extract.ncl"
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

;========================================================
; User options
;========================================================

; unsmoothed ET traj file
etfile = namelist_extract( nlfile, "etfileori", "string")

; type string (used for defining output), generally a model shortname
typestr = namelist_extract( nlfile, "type", "string")

; running average time period (hrs)
hourmean = 24.

; frequency of output data (hrs)
datafreq = namelist_extract( nlfile, "hrintvl", "float")

; average file
catfile = namelist_extract( nlfile, "etfileavg", "string")

;**********************************************************************************

;========================================================
; Read in ET storm trajectories
;========================================================

  etdata_tmp = asciiread(etfile,-1,"string")  ; read in trajectory data as strings
  dims = dimsizes(etdata_tmp)   ; get dimensions
  nrows = dims(0)               ; get number of rows
  ncol = 13                     ; assign number of columns
  delete([/etfile,dims/])

; create new array that will hold all data except rows containing start information
; add column to contain storm number
  etdata_all = new((/nrows,ncol+1/),float)
  
  storm_all = new(nrows,integer)
  
; loop over lines of trajectory data, convert string to float
  n2 = 0       ; counter data without header lines
  num_tot = 0  ; total number of storms
  do n1 = 0,nrows-1  ; loop over all data with header lines
    first = str_get_cols(etdata_tmp(n1),0,0)   ; obtain first character of row
    if (first .ne. "s") then              ; if not 's', this is not a header
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"    "," ")  ; remove extra white space
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"   "," ")   ; remove extra white space
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"  "," ")    ; remove extra white space
      etdata_tmp(n1) = str_strip(etdata_tmp(n1))               ; remove extra white space
      etdata_all(n2,1::) = tofloat(str_split_csv(etdata_tmp(n1)," ",0)) ; convert to float
      etdata_all(n2,0) = s_num   ; assign first dimension as storm number
      n2 = n2 + 1
    else  ; if 's' in the header
      s_num  = stringtoint(str_get_cols(etdata_tmp(n1),41,44))
      storm_all(num_tot) = s_num
      num_tot = num_tot + 1
      continue
    end if
  end do
  
; remove extra lines of missing data
  etdata = etdata_all(0:n2-1,:)
  storm = storm_all(0:num_tot-1)

  delete([/etdata_tmp,etdata_all,storm_all,n1,n2,first,nrows,ncol/])

;========================================================
; Get data and dates
;========================================================

; read in data
  snum = etdata(:,0)
  lon  = etdata(:,1)
  lat  = etdata(:,2)
  pres = etdata(:,3)
  V    = where(etdata(:,4) .eq. -999, etdata@_FillValue,etdata(:,4))
  dist = where(etdata(:,5) .eq. -999, etdata@_FillValue,etdata(:,5))
  ang  = where(etdata(:,6) .eq. -999, etdata@_FillValue,etdata(:,6))
  B    = where(etdata(:,7) .eq. -999, etdata@_FillValue,etdata(:,7))
  Vlt  = where(etdata(:,8) .eq. -999, etdata@_FillValue,etdata(:,8))
  Vut  = where(etdata(:,9) .eq. -999, etdata@_FillValue,etdata(:,9))
  year = etdata(:,10)
  mon  = etdata(:,11)
  day  = etdata(:,12)
  hour = etdata(:,13)
  
  delete([/etdata/])

;========================================================
; Loop over each storm
;========================================================

do i = 0,num_tot-1
  itc = storm(i)   ; storm number
  print("Storm Number: "+itc)

; select range for individual storm
  etsub = ind(snum .eq. itc)
  
; ignore first data point (missing values mess up running avg)
  lon_sub  = lon(etsub(1::))
  lat_sub  = lat(etsub(1::))
  pres_sub = pres(etsub(1::))
  V_sub    = V(etsub(1::))
  dist_sub = dist(etsub(1::))
  ang_sub  = ang(etsub(1::))
  B_sub    = B(etsub(1::))
  Vlt_sub  = Vlt(etsub(1::))
  Vut_sub  = Vut(etsub(1::))
  year_sub = year(etsub(1::))
  mon_sub  = mon(etsub(1::))
  day_sub  = day(etsub(1::))
  hour_sub = hour(etsub(1::))
  
  delete(etsub)

;========================================================
; running avg
;========================================================

; perform running avg over hourmean hours
  nave = toint(hourmean/datafreq) + 1
  opt = 1    ; reflective (symmetric) conditions
  
  B_avg = runave_Wrap(B_sub,nave,opt)
  Vlt_avg = runave_Wrap(Vlt_sub,nave,opt)
  Vut_avg = runave_Wrap(Vut_sub,nave,opt)
  
  dur = dimsizes(B_avg)
  
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
