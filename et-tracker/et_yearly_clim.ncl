;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "../functions/mask_tc.ncl"
;========================================================
; reads in TC and ET trajectories
; generate ET climatology text files
; D. Thatcher 2014/07/23
; v0.1/cleanup C. Zarzycki 2017/08/29
;========================================================
begin
print("et_yearly_clim.ncl")
print(" ")

; choose data set to use
  type = "era"   ; uni, var, or rea
  year_min = 1984
  year_max = 2014

; et phase space criteria  
  Blim = 15.   ; limit for B defining transition (typically 10)
  Vlim = 0.   ; limit for Vlt defining transition (typically 0)

  been_tc_thres = 6   ; how many timesteps does it need to have been TC before counting
  numconsec = 3       ; number of consecutive B/VLT flips to consider ET final

; TC trajectory file
  tfile = "./tc-tracking/new.txt"

; ETC trajectory file
  etfile = "./text_files/traj_et_era_avg"


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
; Read in TC-only trajectories from traj file
;========================================================
 
 print("Loading TC trajectories")
 print(" ")
  
; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
  nstorms_tot = stringtoint(systemfunc("grep 'start' "+tfile+" | wc -l"))

  tdata_tmp = asciiread(tfile,-1,"string")  ; read in trajectory data as strings
  dims = dimsizes(tdata_tmp)    ; get dimensions
  nrows = dims(0)               ; get number of rows
  ncol = 8                      ; assign number of columns
  delete([/tfile,dims/])

; create new array that will hold all data except rows containing start information
; add column to contain storm number
  tdata = new((/nrows-nstorms_tot,ncol+1/),float)
  
; create new array to hold all header information for each storm
  s_year  = new(nstorms_tot,integer)
  s_month = new(nstorms_tot,integer)
  s_day   = new(nstorms_tot,integer)
  s_num   = new(nstorms_tot,integer)

; loop over lines of trajectory data, convert string to float
  nstorm = 0   ; assign storm number
  n2 = 0       ; counter data without header lines
  do n1 = 0,nrows-1  ; loop over all data with header lines
    first = str_get_cols(tdata_tmp(n1),0,0)   ; obtain first character of row
    if (first .ne. "s") then              ; if not 's', this is not a header
      tdata_tmp(n1) = str_sub_str(tdata_tmp(n1),"    "," ")   ; remove extra white space
      tdata_tmp(n1) = str_sub_str(tdata_tmp(n1),"   "," ")    ; remove extra white space
      tdata_tmp(n1) = str_sub_str(tdata_tmp(n1),"  "," ")     ; remove extra white space
      tdata_tmp(n1) = str_strip(tdata_tmp(n1))                ; remove extra white space
      tdata(n2,1::) = tofloat(str_split_csv(tdata_tmp(n1)," ",0))  ; convert to float
      tdata(n2,0) = s_num(nstorm-1)   ; assign first dimension as storm number
      n2 = n2 + 1
    else  ; if 's' in the header
      s_year(nstorm)  = stringtoint(str_get_field(tdata_tmp(n1),3," "))
      s_month(nstorm) = stringtoint(str_get_field(tdata_tmp(n1),4," "))
      s_day(nstorm)   = stringtoint(str_get_field(tdata_tmp(n1),5," "))
      s_num(nstorm)  =  stringtoint(str_get_field(tdata_tmp(n1),7," "))
      nstorm = nstorm + 1
      continue
    end if
  end do

  delete([/tdata_tmp,n1,n2,nstorms_tot,nstorm,first,nrows,ncol/])

; data now contains all trajectory information with the following columns
; 0 = storm number
; 1 = longitude
; 2 = latitude
; 3 = wind speed
; 4 = sea level pressure
; 5 = year
; 6 = month
; 7 = day
; 8 = hour

;========================================================
; Read in ET storm trajectories
;========================================================

print("Loading ET trajectories")
print(" ")

  etdata_tmp = asciiread(etfile,-1,"string")  ; read in trajectory data as strings
  dims = dimsizes(etdata_tmp)   ; get dimensions
  nrows = dims(0)               ; get number of rows
  ncol = 13                     ; assign number of columns
  delete([/etfile,dims/])

; create new array that will hold all data except rows containing start information
; add column to contain storm number
  etdata = new((/nrows,ncol+1/),float)
  
; loop over lines of trajectory data, convert string to float
  n2 = 0       ; counter data without header lines
  do n1 = 0,nrows-1  ; loop over all data with header lines
    first = str_get_cols(etdata_tmp(n1),0,0)   ; obtain first character of row
    if (first .ne. "s") then              ; if not 's', this is not a header
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"    "," ")  ; remove extra white space
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"   "," ")   ; remove extra white space
      etdata_tmp(n1) = str_sub_str(etdata_tmp(n1),"  "," ")    ; remove extra white space
      etdata_tmp(n1) = str_strip(etdata_tmp(n1))               ; remove extra white space
      etdata(n2,1::) = tofloat(str_split_csv(etdata_tmp(n1)," ",0))  ; convert to float
      etdata(n2,0) = snum   ; assign first dimension as storm number
      n2 = n2 + 1
    else  ; if 's' in the header
      snum  = stringtoint(str_get_cols(etdata_tmp(n1),41,44))
      continue
    end if
  end do

; define variables using etdata
  etnum = etdata(:,0)
  psl   = etdata(:,3)
  wind  = etdata(:,4)
;  dist   = etdata(:,5)
;  ang  = etdata(:,6)
  B     = etdata(:,7)
  Vlt   = etdata(:,8)
  year  = etdata(:,10)
  mon   = etdata(:,11)
  day   = etdata(:,12)
  hr    = etdata(:,13)

  delete([/etdata_tmp,n1,n2,first,nrows,ncol,etdata,snum/])

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
  print("Storm # "+itc)

; select trajectory information for individual storm
  tsub_ind = ind(tdata(:,0) .eq. itc)
  tsub = tdata(tsub_ind,:)
  
  dims = dimsizes(tsub)
  nrows = dims(0)

  delete([/tsub_ind,dims/])

;========================================================
; Test for correct basin
;========================================================

; select appropriate basin based on trajectory lat/lon
; analyze only north atlantic basin storms, skip all others
  
  basin = 1
  maskspacing=1.0
  tcMask = mask_tc(0,-999,-999,maskspacing,360)

  orilat=tsub(0,2)
  orilon=tsub(0,1)
  if (orilon .gt. 360.0-maskspacing) then
    orilon=360.0-maskspacing
  end if
  
  if (tcMask({orilat},{orilon}) .ne. basin) then
	delete([/tsub,nrows,basin,maskspacing,tcMask,orilat,orilon/])
	print("  Wrong basin")
	continue
  else if (tsub(nrows-1,2) .lt. 20) then
    n_atl = n_atl + 1
    n_tc = n_tc + 1
    m_atl = m_atl + 1
    m_tc = m_tc + 1
    print("Storm Number: "+itc+"  TC  "+s_year(jtc)+" "+s_month(jtc)+" "+s_day(jtc))
	delete([/tsub,nrows,basin,maskspacing,tcMask/])
	continue
  end if
  end if

  delete([/tsub,basin,maskspacing,tcMask/])

  n_atl = n_atl + 1
  m_atl = m_atl + 1

;========================================================
; Check if ET occurs
;========================================================

; for storm in Atlantic basin, get subset of data 
  et_ind = ind(etnum .eq. itc)
  B_tmp = B(et_ind)
  Vlt_tmp = Vlt(et_ind)
  psl_tmp = psl(et_ind)
  
  year_tmp = year(et_ind)
  mon_tmp = mon(et_ind)
  day_tmp = day(et_ind)
  hr_tmp = hr(et_ind)
  dur = dimsizes(et_ind)

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

  do i = 0,dur-1

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

  delete([/wc,cc,et,no,et_ind,B_tmp,Vlt_tmp,psl_tmp,year_tmp,mon_tmp,day_tmp,hr_tmp,dur/])

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

end
