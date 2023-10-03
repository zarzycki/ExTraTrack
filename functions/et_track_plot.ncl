undef("et_track_plot")
procedure et_track_plot(
  type[1]    : string,   ; type of data
  gifsdir[1] : string,   ; base directory to store images
  itc[1]     : integer,  ; storm number index
  i[1]       : integer,  ; time index
  lat        : numeric,  ; latitudes
  lon        : numeric,  ; longitudes
  date[1]    : numeric,  ; date (YYYYMMDD)
  minP[1]    : numeric,  ; minimum pressure (Pa)
  minPlat    : numeric,  ; latitude of minimum pressure
  minPlon    : numeric,  ; longitude of minimum pressure
  data1      : numeric,  ; data for top panel
  data2      : numeric,  ; data for bottom panel
  debug      : logical   ; if true, optional debugging info
  )
local plot,wks,res,txid0,txid1,latmin,latmax, \
      lon2,lonmin,lonmax,minPlon_ind,minPlon2,data1_flip,data2_flip
; plots the sea level pressure and cyclone center
;
; D.R. Thatcher 15/07/2014
begin

  system("mkdir -p "+gifsdir+"/gifs")

; find mins and maxes for plotting
  latmin = min(lat)
  latmax = max(lat)
  lonmin = min(lon)
  lonmax = max(lon)

; plot storm
  plot = new(2,graphic)
  outfilename=gifsdir+"/gifs/track_plot_"+type+"_"+sprinti("%0.4i",itc)+ \
                      "_"+sprinti("%0.4i",i)
  if (debug) then
    print("Writing: "+outfilename)
  end if
  wks = gsn_open_wks("pdf",outfilename)
  gsn_define_colormap(wks,"rainbow")

  res = True
  res@gsnDraw = False
  res@gsnFrame = False

; fixed contour range
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = 980.0
  res@cnMaxLevelValF = 1030.0
  res@cnLevelSpacingF = 5.0

; add map
  res@gsnAddCyclic = False
  res@mpMinLonF = lonmin
  res@mpMaxLonF = lonmax
  res@mpMinLatF = latmin
  res@mpMaxLatF = latmax
  res@mpCenterLonF = lonmin+(lonmax-lonmin)/2

; include date and minimum pressure value
  res@gsnLeftString = date
  res@gsnRightString = "Min. pressure = "+minP+" Pa"

  plot(0) = gsn_csm_contour_map(wks,data1,res)
  plot(1) = gsn_csm_contour_map(wks,data2,res)

; overlay x at location of minimum pressure
  txid0 = gsn_add_text(wks,plot(0),"x",minPlon,minPlat,False)
  ;txid1 = gsn_add_text(wks,plot(1),"x",minPlon,minPlat,False)

  gsn_panel(wks,plot,(/2,1/),False)

  ;delete([/plot,wks,res,txid0,txid1/])
  delete([/plot,wks,res,txid0/])

end
