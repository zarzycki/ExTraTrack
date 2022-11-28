;========================================================
load "../functions/mask_tc.ncl"
load "../functions/namelist_extract.ncl"
load "../functions/getTrajectoriesGeneric.ncl"
load "../functions/getStormIndices.ncl"
;========================================================
; reads in TC and ET trajectories
; generate ET climatology text files
; D. Thatcher 2014/07/23
; v0.1/cleanup C. Zarzycki 2017/08/29
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
basin = namelist_extract( nlfile, "basin", "integer")
etfile = namelist_extract( nlfile, "etfileavg", "string")
tfile = namelist_extract( nlfile, "tfile", "string")

nldefaults="./defaults/defaults.nl"
maskspacing = namelist_extract( nldefaults, "maskspacing", "float")
; limit for B defining transition (typically 10)
Blim = namelist_extract( nldefaults, "Blim", "float")
; limit for Vlt defining transition (typically 0)
Vlim = namelist_extract( nldefaults, "Vlim", "float")
; how many timesteps does it need to have been TC before counting
been_tc_thres = namelist_extract( nldefaults, "been_tc_thres", "integer")
 ; number of consecutive B/VLT flips to consider ET final
numconsec = namelist_extract( nldefaults, "numconsec", "integer")

; set file names
outfile    = "storms_"+type
monfile    = "mon_"+type
monthfile  = "monthly_"+type
yearfile   = "yearly_"+type
etpathfile = "etpath_yearly_"+type
lifefile   = "life_"+type
etdetailsfile= "etdetails_"+type

;========================================================
; Make subdir for climatology files
;========================================================

 climodir="./climatology_files/"
 system("mkdir -p "+climodir)

;========================================================
; Generate tc mask for masking basins
;========================================================

  tcMask = mask_tc(0,-999,-999,maskspacing,360)

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

delete([/tdata,nstorms_tot,headerDelimStr,maxLines,ncols/])


  
;========================================================
; Count storm types
;========================================================

; loop over each year  
do iyear = year_min,year_max

  n_atl = 0  ; number of atlantic basin storms
  n_tc = 0   ; number of storms that remained tropical
  n_wc = 0   ; number of storms that dissipated as warm cores
  n_cc = 0   ; number of storms that dissipated as cold cores
  n_et = 0   ; number of storms that complete transition
  n_no = 0   ; number of storms that complete partial transition
  n_miss = 0 ; number of missing storms

  et_path1=0
  et_path2=0
  et_path3=0

; remove file for storm numbers of ET
  system("rm -f "+climodir+"/"+outfile+"_"+iyear+".txt")

; loop over each month
do imon = 1,12

; prepare file for output
  mon_str = sprinti("%0.2i",imon)
  system("rm -f ./climatology_files/"+monfile+"_"+iyear+"_"+mon_str+".txt")

  m_atl = 0  ; number of atlantic basin storms
  m_tc = 0   ; number of storms that remained tropical
  m_wc = 0   ; number of storms that dissipated as warm cores
  m_cc = 0   ; number of storms that dissipated as cold cores
  m_et = 0   ; number of storms that complete transition
  m_no = 0   ; number of storms that complete partial transition

;========================================================
; Select storms within desired time period
;========================================================
  
; find indices where data is between given year/month min and max
  ind_storm = ind(s_year .eq. iyear .and. s_month .eq. imon)
  ntc = dimsizes(ind_storm)
  
;========================================================
; Read in data
;========================================================

if (ismissing(ind_storm(0))) then

; write output to text file and remove text file if it already exists
  write_table(climodir+"/"+monfile+"_"+iyear+"_"+mon_str+".txt",\
        "w",[/iyear,imon,m_atl,m_tc,m_wc,m_cc,m_et,m_no/], \
        "%04i  %02i   %02i   %02i   %02i   %02i   %02i   %02i")
        
  delete([/ind_storm/])
  continue
  
else

s_num_sub = s_num(ind_storm)

; loop over each storm
do loop_ind = 0,ntc-1
  itc = s_num_sub(loop_ind)
  jtc = ind(s_num .eq. itc)
  
  print(" ")
  print("Storm # "+itc+" "+jtc)
  
;========================================================
; Test for correct basin
;========================================================

; select appropriate basin based on trajectory lat/lon
; analyze only north atlantic basin storms, skip all others

  orilat=t_lat(jtc,0)
  orilon=t_lon(jtc,0)
  if (orilon .gt. 360.0-maskspacing) then
    orilon=360.0-maskspacing
  end if
  
  print("-- storm originating at lat: "+orilat+"  lon: "+orilon)
    
  if (basin .gt. 0 .and. tcMask({orilat},{orilon}) .ne. basin) then
    delete([/orilat,orilon/])
    print("  Wrong basin")
    continue
  else if (abs(t_lat(jtc,s_dur(jtc)-1)) .lt. 20) then
    n_atl = n_atl + 1
    n_tc = n_tc + 1
    m_atl = m_atl + 1
    m_tc = m_tc + 1
    print("Storm Number: "+itc+"  TC  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
	  continue
  end if
  end if

  n_atl = n_atl + 1
  m_atl = m_atl + 1

;========================================================
; Check if ET occurs
;========================================================

  dur = s_dur(jtc)
  
  B_tmp = B(jtc,0:dur-1)
  Vlt_tmp = Vlt(jtc,0:dur-1)
  psl_tmp = psl(jtc,0:dur-1)
  
  year_tmp = t_yyyy(jtc,0:dur-1)
  mon_tmp = t_mm(jtc,0:dur-1)
  day_tmp = t_dd(jtc,0:dur-1)
  hr_tmp = t_hh(jtc,0:dur-1)

; predefine logicals and et clock
  wc = False   ; warm core dissipation
  cc = False   ; cold core dissipation
  et = False   ; et occurs
  no = False   ; storm isn't a TC long enough in CPS to count
  tc = False   ; storm ends as TC

  et_occurs = False

  et_b=False
  et_v=False

  been_tc = False        ; start tc boolean false
  been_tc_count = 0      ; start count for TC check (should be zero)

  ; start at i = 1 because the first data point has missing values
  do i = 1,dur-1

    ; First, check that storm has been TC for at least X timesteps
    if ( .not. been_tc .and. B_tmp(i) .lt. Blim .and. Vlt_tmp(i) .gt. Vlim) then
      been_tc_count = been_tc_count + 1
      if (been_tc_count .ge. been_tc_thres) then
	      been_tc = True
      end if
    end if

    if (been_tc) then
      if (.not. et_b .and. B_tmp(i) .ge. Blim) then 
        et_b = True
        ; check next 4 things
        jjend=min((/dur-1,(i+numconsec)-1/))
        do jj = i,jjend
          print(" "+B_tmp(jj))
          if (B_tmp(jj) .lt. Blim)
            et_b = False
            print("False start on B")
            break
          end if
        end do
        if (et_b) then
          etb = i
          print("Transition (B) "+year_tmp(i)+" "+mon_tmp(i)+" "+day_tmp(i)+" "+hr_tmp(i))
        end if
      end if

      if (.not. et_v .and. Vlt_tmp(i) .le. Vlim) then
        et_v = True
        ; check next 4 things
        jjend=min((/dur-1,(i+numconsec)-1/))
        do jj = i,jjend
          print(" "+Vlt_tmp(jj))
          if (Vlt_tmp(jj) .gt. Vlim)
            et_v = False
            print("False start on VT")
            break
          end if
        end do
        if (et_v) then
          etv = i
          print("Transition (VT) "+year_tmp(i)+" "+mon_tmp(i)+" "+day_tmp(i)+" "+hr_tmp(i))
        end if
      end if

      ;if (.not. et_v .and. Vlt_tmp(i) .le. Vlim) then
      ;  if ((i .eq. dur-1) .or. (Vlt_tmp(i+1) .le. Vlim)) then
      ;    et_v = True
      ;    etv = i
      ;    print("Transition (VT) "+year_tmp(i)+" "+mon_tmp(i)+" "+day_tmp(i)+" "+hr_tmp(i))
      ;  else
      ;    et_v = False
      ;  end if
      ;end if

      ; check for multiple transitions 
      if (et_b .and. et_v .and. \
        B_tmp(i) .lt. Blim .and. Vlt_tmp(i) .gt. Vlim) then
        print("Warning: could be multiple transitions!")
      end if

    end if
  end do  ; loop over time steps

  whichPath=0  ; flag for which path storm takes to ET
  if (et_b .and. et_v) then
    et_occurs=True
    et = True
    if (etb .lt. etv) then
      ets=etb
      ete=etv
      et_path1=et_path1+1
      whichPath=1
      print("*** ET via traditional pathway")
    else if (etb .gt. etv) then
      ets=etv
      ete=etb
      et_path2=et_path2+1
      whichPath=2
      print("*** non-traditional ET")
    else
      ets=etb
      ete=etv
      et_path3=et_path3+1
      whichPath=3
      print("*** et direct to cold core")
    end if
    end if
  else if (et_b .and. .not. et_v) then
    wc = True
  else if (.not. et_b .and. et_v) then
    cc = True
  else
    if (been_tc) then
      tc = True
    else
      no = True
    end if
  end if
  end if
  end if

  if (et_occurs) then
    print("   Transition completed")
    ; write times to output file
    system("rm -f "+climodir+"/"+lifefile+"_"+sprinti("%0.4i",itc)+".txt")
      write_table(climodir+"/"+lifefile+"_"+sprinti("%0.4i",itc)+".txt","w", \
      [/"storm",itc,dur,year_tmp(0),mon_tmp(0),day_tmp(0),hr_tmp(0)/], \
      "%s    %04i   %03i    %04.0f  %02.0f  %02.0f  %02.0f")
        
       write_table(climodir+"/"+lifefile+"_"+sprinti("%0.4i",itc)+".txt","a", \
       [/itc,ete-ets,year_tmp(ets),mon_tmp(ets),day_tmp(ets),hr_tmp(ets), \
                   year_tmp(ete),mon_tmp(ete),day_tmp(ete),hr_tmp(ete)/], \
        "         %04i   %03i    %04.0f  %02.0f  %02.0f  %02.0f"+\
                       "    %04.0f  %02.0f  %02.0f  %02.0f")

    ; Print each storm's ETC details (colin added)
    system("rm -f "+climodir+"/"+etdetailsfile+"_"+sprinti("%0.4i",itc)+".txt")
      write_table(climodir+"/"+etdetailsfile+"_"+sprinti("%0.4i",itc)+".txt","w", \
       [/itc,ete-ets,whichPath,psl_tmp(ets),psl_tmp(ete),year_tmp(ets),mon_tmp(ets),day_tmp(ets),hr_tmp(ets), \
                   year_tmp(ete),mon_tmp(ete),day_tmp(ete),hr_tmp(ete)/], \
        " %04i   %03i  %01i  %7.2f %7.2f   %04.0f  %02.0f  %02.0f  %02.0f"+\
                       "    %04.0f  %02.0f  %02.0f  %02.0f")
  
  delete(whichPath)
  end if


; count number of each type
  if (et) then
    n_et = n_et + 1
    m_et = m_et + 1
    print("Storm Number: "+itc+"  ET  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
    write_table(climodir+"/"+outfile+"_"+iyear+".txt", \
                "a",[/itc/],"%04i")
  else if (wc) then
    n_wc = n_wc + 1
    m_wc = m_wc + 1
    print("Storm Number: "+itc+"  WC  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
  else if (cc) then
    n_cc = n_cc + 1
    m_cc = m_cc + 1
    print("Storm Number: "+itc+"  CC  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
  else if (no) then
    n_no = n_no + 1
    m_no = m_no + 1
    print("Storm Number: "+itc+"  NO  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
  else if (.not. wc .and. .not. cc .and. .not. et) then
    n_tc = n_tc + 1
    m_tc = m_tc + 1
    print("Storm Number: "+itc+"  TC  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
  end if
  end if
  end if
  end if
  end if

  delete([/wc,cc,et,no,B_tmp,Vlt_tmp,psl_tmp,year_tmp,mon_tmp,day_tmp,hr_tmp,dur/])

end do   ; storm loop

  delete([/s_num_sub/])

end if   ; end if statement for no storms in a given month

;========================================================
; Write output
;========================================================

; write output to text file and remove text file if it already exists
  write_table(climodir+"/"+monfile+"_"+iyear+"_"+mon_str+".txt",\
        "w",[/iyear,imon,m_atl,m_tc,m_wc,m_cc,m_et,m_no/], \
        "%04i  %02i   %02i   %02i   %02i   %02i   %02i   %02i")

  delete([/ind_storm/])

end do   ; month loop

; concatenate monthly output into a single file for each year
  system("cat "+climodir+"/"+monfile+"_"+iyear+\
     "*.txt > "+climodir+"/"+monthfile+"_"+iyear+".txt")
      
; remove excess files
  system("rm -f "+climodir+"/"+monfile+"_"+iyear+"*.txt")

; print output to screen
  print(" ")
  print("Year = "+iyear)
  print("   Atlantic storms: "+n_atl)
  print("    TC dissipators: "+n_tc)
  print("    WC dissipators: "+n_wc)
  print("    CC dissipators: "+n_cc)
  print("    ET transitions: "+n_et)
  print("    NO transitions: "+n_no)
  if (n_miss .gt. 0) then
    print("WARNING, MISSING ET STATS: "+n_miss)
  end if

  print(" ")

; write output to text file and remove text file if it already exists
  system("rm -f "+climodir+"/"+yearfile+"_"+iyear+".txt")
  write_table(climodir+"/"+yearfile+"_"+iyear+".txt",\
        "w",[/iyear,n_atl,n_tc,n_wc,n_cc,n_et,n_no/], \
        "%04i   %02i   %02i   %02i   %02i   %02i   %02i")
              
  delete([/n_et,n_wc,n_cc,n_tc,n_no/])

; write output to text file and remove text file if it already exists
  system("rm -f "+climodir+"/"+etpathfile+"_"+iyear+".txt")
  write_table(climodir+"/"+etpathfile+"_"+iyear+".txt",\
        "w",[/iyear,et_path1,et_path2,et_path3/], \
        "%04i   %02i   %02i   %02i")

  delete([/et_path1,et_path2,et_path3/])

end do   ; year loop

; concatenate monthly files for all years and remove old files
  system("cat "+climodir+"/"+monthfile+"_"+"{"+year_min+".."+year_max+"}"+ \
    ".txt > "+climodir+"/"+monthfile+"_"+year_min+"_"+year_max+".txt")
  system("rm -f ./climatology_files/"+monthfile+"_"+"{"+year_min+".."+year_max+"}"+".txt")

; concatenate yearly files for all years and remove old files
  system("cat "+climodir+"/"+yearfile+"_"+"{"+year_min+".."+year_max+"}"+ \
    ".txt > "+climodir+"/"+yearfile+"_"+year_min+"_"+year_max+".txt")
  system("rm -f "+climodir+"/"+yearfile+"_"+"{"+year_min+".."+year_max+"}"+".txt")

; concatenate etpathfile files for all years and remove old files
  system("cat "+climodir+"/"+etpathfile+"_"+"{"+year_min+".."+year_max+"}"+ \
    ".txt > "+climodir+"/"+etpathfile+"_"+year_min+"_"+year_max+".txt")
  system("rm -f "+climodir+"/"+etpathfile+"_"+"{"+year_min+".."+year_max+"}"+".txt")
     
; concatenate storm files for all years and remove old files
  system("cat "+climodir+"/"+outfile+"_"+"{"+year_min+".."+year_max+"}"+ \
    ".txt > "+climodir+"/"+outfile+".txt")
  system("rm -f "+climodir+"/"+outfile+"_"+"{"+year_min+".."+year_max+"}"+".txt")

; concatenate lifecycle files for all years and remove old files     
  system("cat "+climodir+"/"+lifefile+"_????.txt > "+\
             climodir+"/"+lifefile+".txt")
  system("rm -f "+climodir+"/"+lifefile+"_????.txt")

; concatenate et_detail files for all years and remove old files     
  system("cat "+climodir+"/"+etdetailsfile+"_????.txt > "+\
             climodir+"/"+etdetailsfile+".txt")
  system("rm -f "+climodir+"/"+etdetailsfile+"_????.txt")

  ; final cleanup
  delete([/basin,maskspacing,tcMask/])

end
