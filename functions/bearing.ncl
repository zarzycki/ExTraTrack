undef("bearing")
function bearing(
  lat1[1] : numeric,   ; lat1 in degrees
  lon1[1] : numeric,   ; lon1 in degrees
  lat2[1] : numeric,   ; lat2 in degrees
  lon2[1] : numeric    ; lon2 in degrees
  )
local pi, d2r, r2d, lat1r, lon1r, lat2r, lon2r, ang_tmp, lat1_tmp, lon1_tmp
; calculates the bearing (angle of movement over the sphere)
; between two points on a sphere
;
; output angle in degrees
;
; D.R. Thatcher 15/07/2014
begin

; constants
  pi = 4.*atan(1.)
  d2r = pi/180.
  r2d = (1/d2r)

; convert lat/lon points from degrees to radians
  lat1r = lat1*d2r
  lon1r = lon1*d2r
  lat2r = lat2*d2r
  lon2r = lon2*d2r
  
; compute angle of motion based on two lat/lon points
  ang_tmp = r2d*atan2(sin((lon2r-lon1r))*cos(lat2r), \
                  cos(lat1r)*sin(lat2r) - sin(lat1r)*cos(lat2r)*cos(lat2r-lat1r))

; convert angles from (-180 to +180) to (0 to 360)
  ang = mod(ang_tmp + 360.,360.)

  delete([/pi, d2r, r2d, lat1r, lon1r, lat2r, lon2r, ang_tmp/])
               
  return(ang)
  
end