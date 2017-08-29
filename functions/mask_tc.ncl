function mask_tc(gridchoice:integer,latin[*]:numeric,lonin[*]:numeric,gridsize[*]:numeric,lonstart:integer)

; Gridchoice -- specify how you want to generate the TC mask
; 0 means you specify a regular lat/lon spacing (gridsize) and a lonstart
; 1 means you pass in lat/lon (latin/lonin) arrays
; latin - 
; lonin - 
; gridsize - grid spacing
; lon

begin

if (gridchoice .eq. 0) then
  latS = -90
  latN =  90
  lonW =   0
  lonE = 360

  dlat =  gridsize
  dlon =  gridsize

  nlat = toint((latN-latS)/dlat) + 1
  mlon = toint((lonE-lonW)/dlon)

  lat  = fspan(latS, latN, nlat)
  lon  = fspan(lonW, lonE-dlon, mlon)
  lat@units = "degrees_north"
  lon@units = "degrees_east"
else if (gridchoice .eq. 1) then
  lat = latin
  lon = lonin
  nlat = dimsizes(lat)
  mlon = dimsizes(lon)
  lat@units = "degrees_north"
  lon@units = "degrees_east"
else
  print("Invalid gridchoice value")
  exit
end if
end if
  
  maskArr = new((/nlat,mlon/),integer)
  
  maskArr!0 = "lat"
  maskArr!1 = "lon"
  maskArr&lat = lat
  maskArr&lon = lon
  
  ; Coefficients for calculating ATL/EPAC sloped line
  m = -0.58
  b = 00 - m*295
  maxlat = 45.0
 
  ; 1 = ATL
  ; 2 = EPAC
  ; 3 = CPAC
  ; 4 = WPAC
  ; 5 = NIO
  ; 6 = SIO
  ; 7 = SPAC
  ; 8 = SATL
  ; Atlantic
  do i = 0,nlat-1
    do j = 0,mlon-1
      if (lat(i) .ge. 0 .and. lat(i) .le. maxlat .and. lon(j) .gt. 257 .and. lon(j) .le. 355) then
        ;print(lon(j))
        funcval = m*lon(j) + b
      ;print(funcval)
        if (lat(i) .gt. funcval) then
        maskArr(i,j) = 1
        else
        maskArr(i,j) = 2
        end if
      else if (lat(i) .ge. 0 .and. lat(i) .le. maxlat .and. lon(j) .gt. 220 .and. lon(j) .le. 257) then
              maskArr(i,j) = 2
      else if (lat(i) .ge. 0 .and. lat(i) .le. maxlat .and. lon(j) .gt. 180 .and. lon(j) .le. 220) then
              maskArr(i,j) = 3
      else if (lat(i) .ge. 0 .and. lat(i) .le. maxlat .and. lon(j) .gt. 100 .and. lon(j) .le. 180) then
              maskArr(i,j) = 4
      else if (lat(i) .ge. 0 .and. lat(i) .le. maxlat .and. lon(j) .gt. 30 .and. lon(j) .le. 100) then
              maskArr(i,j) = 5
      else if (lat(i) .lt. 0 .and. lat(i) .ge. -maxlat .and. lon(j) .gt. 30 .and. lon(j) .le. 135) then
              maskArr(i,j) = 6
      else if (lat(i) .lt. 0 .and. lat(i) .ge. -maxlat .and. lon(j) .gt. 135 .and. lon(j) .le. 290) then
              maskArr(i,j) = 7   
      
      
      else
        maskArr(i,j) = 0
      end if
      end if
      end if
      end if
      end if
      end if
      end if
    end do
  end do
  
  if (lonstart .ne. 360) then
    maskArr = lonPivot(maskArr,lonstart)
  end if
  
  return(maskArr)

; ;**************************************************
; ; create plot
; ;**************************************************
;   wks = gsn_open_wks("pdf","tc_mask")
;  
; 	res                      = True
; 	res@gsnMaximize          = True
; 	  
; 	res@cnFillOn             = True              ; turn on color
;   res@cnLinesOn            = False             ; turn off contour lines
;   res@cnLineLabelsOn       = False
;  
;   res@gsnAddCyclic = False
;   
;   plot = gsn_csm_contour_map(wks,maskArr,res)

end

