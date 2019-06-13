;========================================================
; Track density
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;========================================================

; plot AMIP trajectories

begin

  ;type = "rea"
  dat = "avg"
  plot_legend = False
  file_type = "pdf"



  itc = toint(itcstr)
  ;yyyy = "2002"
  ;mm1 = "09"
  ;mm2 = "11"
  ;date_start = "2002 10 16 12"
  ;date_end = "2002 10 17 12"


  print(" ")
  print("Storm Number:  " + itc)
  


;========================================================
; Read in surface pressure
;========================================================

;   fnam1 = "/global/homes/c/czarzyck/scratch/et-hyperion/"+\
;           "*."+yyyy+"-{"+mm1+".."+mm2+"}*"
;   files1 = systemfunc("ls "+fnam1)
;   input = addfiles(files1,"r")
;       ListSetType(input,"cat")

;  ps_tmp  = input[:]->PSL
;  ps_in = ps_tmp(:,:,:)
;  ps = ps_in/100.0
;    copy_VarCoords(ps_in,ps)
;  date_tmp = cd_calendar(input[:]->time,-2)
;  date_ps = date_tmp(:)
;  
;  delete([/ps_tmp,date_tmp,fnam1,files1,input,ps_in/])

  
;========================================================
; Read in ET storm trajectories
;========================================================

; trajectory directory
  etdir = "../et-tracker/text_files/"
  etfile = "traj_et_"+type+"_"+dat

; read in trajectory data from ascii file
; get number of storms by counting number of lines with 'start'
  nstorms_tot = stringtoint(systemfunc("grep 'start' "+etdir+etfile+" | wc -l"))
  
  etdata_tmp = asciiread(etdir+etfile,-1,"string")  ; read in trajectory data as strings
  dims = dimsizes(etdata_tmp)   ; get dimensions
  nrows = dims(0)               ; get number of rows
  ncol = 13                     ; assign number of columns
  delete([/etdir,etfile,dims/])

; create new array that will hold all data except rows containing start information
; add column to contain storm number
  etdata = new((/nrows-nstorms_tot,ncol+1/),float)
  
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
      etdata(n2,0) = s_num   ; assign first dimension as storm number
      n2 = n2 + 1
    else  ; if 's' in the header
      s_num  = stringtoint(str_get_cols(etdata_tmp(n1),41,44))
      continue
    end if
  end do

  snum = etdata(:,0)
  lon  = etdata(:,1)
  lat  = etdata(:,2)
  p    = etdata(:,3)
  B    = etdata(:,7)
  Vlt  = etdata(:,8)
  Vut  = etdata(:,9)
  snum = etdata(:,0)

  lon@units="degrees_east"
  
; date
  year  = sprintf("%04.0f",etdata(:,10))
  month = sprintf("%02.0f",etdata(:,11))
  day   = sprintf("%02.0f",etdata(:,12))
  hour  = sprintf("%02.0f",etdata(:,13))
  date  = year+" "+month+" "+day+" "+hour
  date_et = year+month+day
  
  delete([/etdata,year,month,day,s_num/])

;========================================================
; Load list of ET storms
;========================================================

; file with storm index for ET
  storm_dir = "../et-tracker/climatology_files/"
  storm_fil = "storms_"+type+".txt"
  
; list of storm numbers  
  storm_num = asciiread(storm_dir+storm_fil,-1,"integer")
  
; total number of ET
  net = dimsizes(storm_num)

; find correct storm
  sub = ind(snum .eq. itc)
  n_itc = dimsizes(sub)

;========================================================
; draw map
;========================================================

  system("mkdir -p phase_traj")
  
  latS = 0
  latN = 80
  lonW = 0
  lonE = 360
    
  i = n_itc-1

    wks  = gsn_open_wks(file_type,"./phase_traj/trajs_"+type+"_"+sprinti("%0.4i",itc)+"_"+sprinti("%0.4i",i))
 
;========================================================
; draw phase space
;========================================================
 
    p_sub = p(sub(0:i))
    B_sub = B(sub(0:i))
    Vlt_sub = Vlt(sub(0:i))
    Vut_sub = Vut(sub(0:i))
    date_sub = date(sub(0:i))
    date_et_sub = date_et(sub(i))
    hour_sub = hour(sub(i))

	plot = new(2,graphic)

    p_min = 930
    p_max = 1005
    p_int = 5

	res           = True
	res@gsnDraw   = False
	res@gsnFrame  = False

	res@tmXBLabelFontHeightF   = 0.02
	res@tmYLLabelFontHeightF   = 0.02
	res@tiXAxisFontHeightF     = 0.025
	res@tiYAxisFontHeightF     = 0.025
	res@gsnXRefLine            = 0
	res@gsnYRefLine            = 15
	res@trXMinF                = -500
	res@trXMaxF                = 500

;========================================================
; First draw lines for both plots
;========================================================

	res@xyMarkLineMode         ="Lines"
	res@xyMarkerColor          = "Black"
  
	res@trYMinF                = -25
	res@trYMaxF                = 100
	res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
	res@tiYAxisString          = "B [900-600 hPa Storm-Relative Thickness Symmetry]"

	plot(0)=gsn_csm_xy(wks,Vlt_sub,B_sub,res)
  
    res@gsnYRefLine = 0
	res@trYMinF                = -400
	res@trYMaxF                = 600
	res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
	res@tiYAxisString          = "-V~B~T~N~~S~U~N~ [600-300 hPa Thermal Wind]"

	plot(1)=gsn_csm_xy(wks,Vlt_sub,Vut_sub,res)
  
	delete(res@gsnXRefLine)
	delete(res@gsnYRefLine)
    
;======================================================== 
; Overlay markers
;========================================================

  plot_int = 4
  start_plot = 1

; dummy data for markers
  dum0 = new(dimsizes(Vlt_sub(start_plot::plot_int)),graphic)
  dum1 = new(dimsizes(Vlt_sub(start_plot::plot_int)),graphic)
  dum2 = new(dimsizes(Vlt_sub(start_plot::plot_int)),graphic)
  dum3 = new(dimsizes(Vlt_sub(start_plot::plot_int)),graphic)

; get info for color map 
  map_str = "BlAqGrYeOrReVi200"
  levels = ispan(p_min,p_max,p_int)
  black = (/0.0,0.0,0.0,1.0/)

  j2 = start_plot  
  do j = 0,dimsizes(Vlt_sub(start_plot::plot_int))-1
  
; get color for marker fill
    rgba = get_color_rgba(map_str,levels,p_sub(j2))
    rgba(0:2) = 1-rgba(0:2)   ; flip colors

    res2 = True
	res2@gsMarkerIndex   = 16
	res2@gsMarkerSizeF   = 0.012
	res2@gsMarkerColor   = rgba

    dum0(j) = gsn_add_polymarker(wks,plot(0),Vlt_sub(j2),B_sub(j2),res2)
    dum1(j) = gsn_add_polymarker(wks,plot(1),Vlt_sub(j2),Vut_sub(j2),res2)

;========================================================
; Overlay marker outline
;========================================================

	res2@gsMarkerIndex   = 4
	res2@gsMarkerSizeF   = 0.012
	res2@gsMarkerColor   = black

    dum2(j) = gsn_add_polymarker(wks,plot(0),Vlt_sub(j2),B_sub(j2),res2)
    dum3(j) = gsn_add_polymarker(wks,plot(1),Vlt_sub(j2),Vut_sub(j2),res2)
  
    j2 = j2 + plot_int
  
  end do   ; end marker loop
  
  j2 = 0  
  do j = 0,dimsizes(Vlt_sub)-1
  
    if (date_sub(j2) .eq. date_start) then
      res2@gsMarkerSizeF = 0.025
      res2@gsMarkerThicknessF = 2
      res2@gsMarkerIndex = 5      
      dum8 = gsn_add_polymarker(wks,plot(0),Vlt_sub(j2),B_sub(j2),res2)
      dum9 = gsn_add_polymarker(wks,plot(1),Vlt_sub(j2),Vut_sub(j2),res2)

    end if
    
    if (date_sub(j2) .eq. date_end) then
      res2@gsMarkerSizeF = 0.025
      res2@gsMarkerThicknessF = 2
      res2@gsMarkerIndex = 6
      dum10 = gsn_add_polymarker(wks,plot(0),Vlt_sub(j2),B_sub(j2),res2)
      dum11 = gsn_add_polymarker(wks,plot(1),Vlt_sub(j2),Vut_sub(j2),res2)
    end if

    j2 = j2 + 1

  end do   ; end marker loop

    resP = True
    resP@gsnFrame = False
    resP@gsnPanelRight = 0.35
    gsn_panel(wks,plot,(/2,1/),resP)
    
    delete([/res,resP/])
    delete([/j,j2,Vlt_sub,B_sub,Vut_sub,res2,dum0,dum1,dum2,dum3/])

;========================================================
; draw trajectory
;========================================================

    ;ps_ind = ind(date_ps .eq. date_et_sub)
   ; 
   ; if (hour_sub .eq. "00") then
   ;   ps_sub = ps(ps_ind(0),:,:)
   ; else if (hour_sub .eq. "06") then
   ;   ps_sub = ps(ps_ind(1),:,:)
   ; else if (hour_sub .eq. "12") then
   ;   ps_sub = ps(ps_ind(2),:,:)
   ; else if (hour_sub .eq. "18") then
   ;   ps_sub = ps(ps_ind(3),:,:)
  ;  end if
  ;  end if
  ;  end if
  ;  end if

    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

	res                       = True     ; plot mods desired
	res@gsnDraw               = False
	res@gsnFrame              = False

	res@mpMinLatF             = latS
	res@mpMaxLatF             = latN
	res@mpMinLonF             = lonW
	res@mpMaxLonF             = lonE
	res@mpCenterLonF          = (lonE+lonW)*0.5

if (type .eq. "rea") then
	res@tiMainString = "CFSR"	
else if (type .eq. "nx2") then
	res@tiMainString = "CAM-55"	
else if (type .eq. "nx4") then
	res@tiMainString = "CAM-28"	
else
	res@tiMainString = ""	
end if
end if
end if

	res@gsnCenterString = date_sub(0)+" hrs to "+date_sub(i)+" hrs"
	
	lineres = True
	lineres@gsLineThicknessF = 5.0

	lon_sub = lon(sub(0:i))
	lat_sub = lat(sub(0:i))

    map = gsn_csm_map_ce(wks,res)
	dum = gsn_add_polyline(wks,map,lon_sub,lat_sub,lineres)
	
;======================================================== 
; Overlay markers
;========================================================

; dummy data for markers
  dum4 = new(dimsizes(lon_sub(start_plot::plot_int)),graphic)
  dum5 = new(dimsizes(lon_sub(start_plot::plot_int)),graphic)

; get info for color map 
  map_str = "BlAqGrYeOrReVi200"
  levels = ispan(p_min,p_max,p_int)
  black = (/0.0,0.0,0.0,1.0/)
   
  j2 = start_plot  
  do j = 0,dimsizes(lon_sub(start_plot::plot_int))-1

; get color for marker fill
    rgba = get_color_rgba(map_str,levels,p_sub(j2))
    rgba(0:2) = 1-rgba(0:2)   ; flip colors
    
    res2 = True
	res2@gsMarkerIndex   = 16
	res2@gsMarkerSizeF   = 0.012
	res2@gsMarkerThicknessF = 1
	res2@gsMarkerColor   = rgba

    dum4(j) = gsn_add_polymarker(wks,map,lon_sub(j2),lat_sub(j2),res2)

;========================================================
; Overlay marker outline
;========================================================

	res2@gsMarkerIndex   = 4
	res2@gsMarkerSizeF   = 0.012
	res2@gsMarkerColor   = black

    dum5(j) = gsn_add_polymarker(wks,map,lon_sub(j2),lat_sub(j2),res2)
  
    j2 = j2 + plot_int
  
  end do   ; end marker loop

  j2 = 0  
  do j = 0,dimsizes(lon_sub)-1
    
;2000  09  28  18    2000  10  01  06
    if (date_sub(j2) .eq. date_start) then
      res2@gsMarkerSizeF = 0.03
      res2@gsMarkerThicknessF = 3
      res2@gsMarkerIndex = 5
      dum6 = gsn_add_polymarker(wks,map,lon_sub(j2),lat_sub(j2),res2)
    end if
    
    if (date_sub(j2) .eq. date_end) then
      res2@gsMarkerSizeF = 0.03
      res2@gsMarkerThicknessF = 3
      res2@gsMarkerIndex = 6
      dum7 = gsn_add_polymarker(wks,map,lon_sub(j2),lat_sub(j2),res2)
    end if
  
    j2 = j2 + 1
  
  end do   ; end marker loop

  y_loc = 0.45
  x_loc = 0.7

  txres = True
  txres@txFontHeightF = 0.015

  res2@gsMarkerSizeF = 0.03
  res2@gsMarkerThicknessF = 3
  res2@gsMarkerIndex = 6
  if (plot_legend) then
    gsn_polymarker_ndc(wks,x_loc,y_loc,res2)
    gsn_text_ndc(wks,"Transition Ends   ",x_loc+0.1,y_loc,txres)
  end if
  res2@gsMarkerIndex = 5
  if (plot_legend) then
    gsn_polymarker_ndc(wks,x_loc,y_loc+0.04,res2)
    gsn_text_ndc(wks,"Transition Begins",x_loc+0.1,y_loc+0.04,txres)
  end if
  labels = levels
  labelcolors = new((/dimsizes(levels),3/),float)
  do n=0,dimsizes(levels)-1
    tmp = get_color_rgba(map_str,levels,levels(n))
    labelcolors(n,:) = 1-tmp(0:2)
  end do


  lbres = True
  lbres@lbPerimOn            = False               ; no label bar box
  lbres@lbOrientation        = "Horizontal"        ; orientation
  lbres@vpWidthF             = 0.6                 ; size
  lbres@vpHeightF            = 0.1
  lbres@lbLabelFontHeightF   = 0.012               ; label font height
  lbres@lbLabelAlignment     = "InteriorEdges"     ; where to label
  lbres@lbMonoFillPattern    = True                ; fill sold
  lbres@lbFillColors         = labelcolors      ; must be RGB triplets
  
  lbres@lbTitleOn = True
  lbres@lbTitlePosition = "Bottom"
  lbres@lbTitleFontHeightF = 0.015
  lbres@lbTitleString = "Minimum Surface Pressure (hPa)"
  
  gsn_labelbar_ndc (wks,dimsizes(levels),labels,0.4,0.25,lbres)

    resP = True
    resP@gsnFrame = False
    resP@gsnPanelLeft = 0.4
    gsn_panel(wks,map,(/1,1/),resP)
    
    frame(wks)
	
	delete([/lon_sub,lat_sub,date_sub,dum,res,resP/])
  
  exit

end
