;========================================================
; Reads 
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "calc_area.ncl"
load "/glade/u/home/zarzycki/ncl/functions/mask_tc.ncl"
;========================================================
begin


;---------------User settings----------------------------

; False, False for BAMS Clivar paper, True/true for JAMES
plot_1deg_too=False
drawRegLines=False

wind_factor = "log"
traj_dir = "/glade/u/home/zarzycki/tcstats/"
traj_filename = "traj_vr"
if (plot_1deg_too) then
  traj_1deg_dir = traj_dir
  traj_1deg_filename = "traj_1deg"
end if
ibdir = "/glade/u/home/zarzycki/work/ibtracs/"
ibfilename="Allstorms.ibtracs_all.v03r04.nc.nc"
out_type = "pdf"
basin = 1
colorArray = (/"blue","red","springgreen4"/)
markerTypes = (/2,5,1/)

threshWind = 10 ; threshold wind in m/s

;--------------------------------------------------------

ms_to_kts = 1.94384449

;--------------------------------------------------------

; Wind speed power relation
if (wind_factor .eq. "power") then
  z_10 = 10.
  z_r = 64.
  alpha = 0.11
  factor = (z_10/z_r)^alpha
  print("Using factor: "+factor)
else if (wind_factor .eq. "log") then
  ; Garret 1992 -- Wind profile formulation
  ; Wieringa 1993 -- roughness coefficient
  z_10 = 10.
  z_r = 64.
  z_0 = 0.0002 ; roughness coeff
  factor = 1 + log(z_10/z_r)/log(z_r/z_0)
  print("Using factor: "+factor)
else if (wind_factor .eq. "none") then
  print("No correction used.")
else
  print("Incorrect wind correction type.")
  exit
end if
end if
end if

print("Getting TC Mask")
  maskspacing=2
  tcMask = mask_tc(0,-999,-999,maskspacing,360)
  tcMaskIB = mask_tc(0,-999,-999,maskspacing,180)

print("Read in TC locations from CAM")
traj_cleanfile = "traj_filt_clean.TMP" ; hardcoded temp file name, feel free to change
print("Removing text headers from trajectories")
system("grep -v '^start' "+ traj_dir+traj_filename +" > "+traj_cleanfile+"")

traj_ncols = numAsciiCol(traj_cleanfile)
traj_data = readAsciiTable(traj_cleanfile,traj_ncols,"float",0)
traj_nrows = dimsizes(traj_data(:,0))    ; calculate # of rows
print("'" + traj_filename + "' has " + traj_nrows + " rows and " + traj_ncols + \
      " columns of data.")
print("Removing temporary file:"+traj_cleanfile)
system("rm -rf "+traj_cleanfile)
        
; wind is column 2, pressure is column 3
  
do i = 0,traj_nrows-1
  cyclon = traj_data(i,0)
  cyclat = traj_data(i,1)
  if (tcMask({cyclat},{cyclon}) .ne. basin)
    ; If not in Atlantic basin, replace with missing values
    traj_data(i,3) = -999
    traj_data(i,2) = -999
  else
  end if
end do

if (threshWind .gt. 0) then
  traj_data(:,3) = where(traj_data(:,3).lt.threshWind, -999, traj_data(:,3))
  traj_data(:,2) = where(traj_data(:,2).lt.threshWind, -999, traj_data(:,2))  
end if

if (plot_1deg_too) then
  traj_cleanfile = "traj_filt_clean.TMP" ; hardcoded temp file name, feel free to change
  print("Removing text headers from trajectories")
  system("grep -v '^start' "+ traj_1deg_dir+traj_1deg_filename +" > "+traj_cleanfile+"")
  
  traj_ncols = numAsciiCol(traj_cleanfile)
  traj_data_1deg = readAsciiTable(traj_cleanfile,traj_ncols,"float",0)
  traj_nrows = dimsizes(traj_data_1deg(:,0))    ; calculate # of rows
  print("'" + traj_1deg_filename + "' has " + traj_nrows + " rows and " + traj_ncols + \
        " columns of data.")
  print("Removing temporary file:"+traj_cleanfile)
  system("rm -rf "+traj_cleanfile)
          
  ; wind is column 2, pressure is column 3
    
  do i = 0,traj_nrows-1
    cyclon = traj_data_1deg(i,0)
    cyclat = traj_data_1deg(i,1)
    if (tcMask({cyclat},{cyclon}) .ne. basin)
      ; If not in Atlantic basin, replace with missing values
      traj_data_1deg(i,3) = -999
      traj_data_1deg(i,2) = -999
    else
    end if
  end do

  if (threshWind .gt. 0) then
    traj_data_1deg(:,3) = where(traj_data_1deg(:,3).lt.threshWind, -999, traj_data_1deg(:,3))
    traj_data_1deg(:,2) = where(traj_data_1deg(:,2).lt.threshWind, -999, traj_data_1deg(:,2))  
  end if

end if
;--------------------------------------------------------

print("Get ibTracs data")

ibfile = addfile(ibdir+ibfilename,"r")

iblat = ibfile->lat_for_mapping*.01
iblon = ibfile->lon_for_mapping*.01
; Quick fix for lons that go outside the tcMask
iblon(:,0) = where(iblon(:,0).gt.180.0-maskspacing,180.0-maskspacing,iblon(:,0))
ibpres = ibfile->source_pres(:,:,0)*.1
ibyear = ibfile->season
ibwind = ibfile->source_wind(:,:,0)*.1/ms_to_kts

ibdimsizes = dimsizes(iblat)
ibnumstorms = ibdimsizes(0)

do i = 0,ibnumstorms-1
  if (ibyear(i) .ge. 1980 .and. ibyear(i) .le. 2002 .and. tcMaskIB({iblat(i,0)},{iblon(i,0)}) .eq. basin) then
  ;if (ibyear(i) .ge. 1981 .and. ibyear(i) .le. 1984) then
  else
    ibpres(i,:) = ibpres@_FillValue
    ibwind(i,:) = ibwind@_FillValue
  end if
end do

onedpres = ndtooned(ibpres)
onedwind = ndtooned(ibwind)

if (threshWind .gt. 0) then
  onedpres = where(onedpres.lt.threshWind, -999, onedpres)
  onedwind = where(onedwind.lt.threshWind, -999, onedwind)  
end if

;--------------------------------------------------------
if drawRegLines then
  print("Calculating regression for best-fit line")
  
  polyn = 3
  quad_obs = lspoly(onedwind, onedpres , 1, polyn)    ; all weights are set to one
  ;print("Calculating regression for cam")
  quad_cam = lspoly(traj_data(:,2)*factor, traj_data(:,3) , 1, polyn)
  ;print(quad_cam) ; Y = c(0) + c(1)*x + c(2)*x^2
  regx = fspan(0,100,750)
  camregy = new(dimsizes(regx),float)
  obsregy = new(dimsizes(regx),float)
  camregy = quad_cam(0) + quad_cam(1)*regx + quad_cam(2)*regx^2
  obsregy = quad_obs(0) + quad_obs(1)*regx + quad_obs(2)*regx^2
  
  if (plot_1deg_too) then
    regx_short = fspan(0,45,750)
    quad_cam_1deg = lspoly(traj_data_1deg(:,2)*factor, traj_data_1deg(:,3) , 1, polyn)
    camregy_1deg = new(dimsizes(regx_short),float)
    camregy_1deg = quad_cam_1deg(0) + quad_cam_1deg(1)*regx_short + quad_cam_1deg(2)*regx_short^2
  end if
end if
  
;-------------------------------------------------------- 

if (plot_1deg_too) then
  plotName = "pres_wind"
else
  plotName = "pres_wind_no_1deg"
end if
wks  = gsn_open_wks(out_type,plotName)            ; open a ps file

res  = True
res@gsnMaximize            = True       ; Maximize plot in frame.
res@gsnDraw    = False                         ; don't draw
res@gsnFrame   = False                         ; don't advance frame

res@tiMainString  = "Pressure-Wind Relationship"	; title
res@tiYAxisString = "Wind speed (m/s)" 	; y axis title
res@tiXAxisString = "Surface pressure (hPa)" 	; y axis title

res@tiYAxisFontHeightF = 0.018
res@tiXAxisFontHeightF = 0.018
res@tiMainFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@trXMinF           = 887.
res@trXMaxF           = 1022.
res@trYMinF           = 0.
res@trYMaxF           = 86.

res@xyMarkLineMode = "Markers"


; Draw ibTracks data
res@xyMarkerColors = colorArray(0)
res@xyMarker      = markerTypes(0)
res@xyMarkerSizeF = 0.005
plot = gsn_csm_xy(wks,onedpres,onedwind,res)

; Draw CAM-SE data
res@xyMarkerColors = colorArray(1)
res@xyMarker      = markerTypes(1)
res@xyMarkerSizeF = 0.004
plot2 = gsn_csm_xy(wks,traj_data(:,3),traj_data(:,2)*factor,res)
overlay(plot,plot2)

if (plot_1deg_too) then
  res@xyMarkerColors = "white"
  res@xyMarker      = markerTypes(2)
  res@xyMarkerSizeF = 0.018 ;0.018
  plot5 = gsn_csm_xy(wks,traj_data_1deg(:,3),traj_data_1deg(:,2)*factor,res)
  overlay(plot,plot5)
  res@xyMarkerColors = colorArray(2)
  res@xyMarker      = markerTypes(2)
  res@xyMarkerSizeF = 0.014 ;0.018
  plot6 = gsn_csm_xy(wks,traj_data_1deg(:,3),traj_data_1deg(:,2)*factor,res)
  overlay(plot,plot6)
end if

; Now we will plot the quad fit lines
if drawRegLines then
  res@xyMarker      = 1
  
  res@xyMarkerSizeF = 0.020
  res@xyMarkerColors = colorArray(0)
  plot3 = gsn_csm_xy(wks,obsregy,regx,res)
  overlay(plot,plot3)
  
  res@xyMarkerColors = colorArray(1)
  plot4 = gsn_csm_xy(wks,camregy,regx,res)
  overlay(plot,plot4)
  
  if (plot_1deg_too) then
    res@xyMarkerColors = colorArray(2)
    plot6 = gsn_csm_xy(wks,camregy_1deg,regx_short,res)
    overlay(plot,plot6)
  end if
end if


lgres                    = True
lgres@lgMonoMarkerSize = False
lgres@lgLabelFontHeightF = .18
lgres@vpWidthF           = 0.24
lgres@vpHeightF          = 0.10
lgres@lgPerimColor       = "black"
lgres@lgPerimThicknessF  = 1.0
lgres@lgMonoMarkerIndex = False
lgres@lgItemType         = "Markers"
if (plot_1deg_too) then
  lgres@lgMarkerColors     = colorArray(::-1)
  lgres@lgMarkerIndexes = markerTypes(::-1)
  lgres@lgMarkerSizes = (/0.028,0.005,0.007/)
  lbid = gsn_create_legend(wks,3,(/"CAM-SE Uni. (1~S~o~N~)","CAM-SE V-R (0.25~S~o~N~)","Obs (IBTrACS)"/),lgres)
else
  lgres@lgMarkerColors     = (/colorArray(1),colorArray(0)/)
  lgres@lgMarkerIndexes = (/markerTypes(1),markerTypes(0)/)
  lgres@lgMarkerSizes = (/0.005,0.007/)
  lbid = gsn_create_legend(wks,2,(/"CAM-SE V-R (0.25~S~o~N~)","Obs (IBTrACS)"/),lgres)
end if

amres = True
amres@amParallelPosF   = 0.3		    ; move legend to the right
amres@amOrthogonalPosF = -.4             ; move the legend down
annoid1 = gsn_add_annotation(plot,lbid,amres)   ; add legend to plot

draw(plot)
frame(wks)

end

