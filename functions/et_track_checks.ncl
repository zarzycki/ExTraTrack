undef("et_track_checks")
function et_track_checks(
  latmax[1]  : numeric,    ; domain latitude boundary
  lonmax[1]  : numeric,    ; domain longitude boundary
  minP[1]    : numeric,    ; minimum pressure (hPa)
  minPlat[1] : numeric,    ; latitude of minimum pressure (degrees)
  minPlon[1] : numeric,    ; longitude of minimum pressure (degrees)
  dist[1]    : numeric,    ; distance traveled (km)
  d_ang[1]   : numeric     ; change in bearing of storm
  )
local vel
; checks if the storm can should still be considered
;
; D.R. Thatcher 15/07/2014
begin

  storm_check = True  ; storm position (True = valid position)

  vel = (dist/(60.*60.*6.))*1000.

; if min pressure above 1020 hPa
; from Hart 2003 (Monthly Weather Review)
  if (minP .ge. 1020.0) then
    storm_check = False
  end if
  
; if storm moves faster than 40 m/s
; from Hart 2003 (Monthly Weather Review)
  if (vel .ge. 40.) then
    storm_check = False
  end if
 
; if the storm goes south of 20 degN or hits northern domain boundary
  if (abs(minPlat) .lt. 20. .or. minPlat .ge. latmax-0.001) then
    storm_check = False
  end if
  
; if the storm hits eastern domain boundary
  if (minPlon .ge. lonmax-0.001) then
    storm_check = False
  end if
  
; if the tracked storm has a sharp change in directions
; from Hart 2003 (Monthly Weather Review)
  if (vel .ge. 10. .and. vel .lt. 15.) then
    if (abs(d_ang) .gt. 135.) then
      storm_check = False
    end if
  else if (vel .ge. 15. .and. vel .lt. 20.) then
    if (abs(d_ang) .gt. 90.) then
      storm_check = False
    end if
  else if (vel .ge. 20. .and. vel .lt. 25.) then
    if (abs(d_ang) .gt. 75.) then
      storm_check = False
    end if
  else if (vel .ge. 25. .and. vel .lt. 30.) then
    if (abs(d_ang) .gt. 60.) then
      storm_check = False
    end if
  else if (vel .ge. 30.) then
    if (abs(d_ang) .gt. 45.) then
      storm_check = False
    end if
  end if
  end if
  end if
  end if
  end if
  
  delete(vel)
  
  return(storm_check)
  
end
