undef("extrap_latlon")
function extrap_latlon(
  lat1[1] : numeric,
  lon1[1] : numeric,
  d[1]    : numeric,
  brng[1] : numeric,
  extrap_latlon[1] : logical
  )
local pi, d2r, r2d, lat1r, lon1r, lat2r, lon2r, brngr, lat2, lon2, R, return_var

begin

  ; Define the Earth's radius in kilometers
  R = 6371.0

  pi = 4.*atan(1.)
  d2r = pi/180.
  r2d = (1/d2r)

  brng = brng

  ; Convert initial latitude, longitude, and bearing to radians
  lat1r = lat1 * d2r
  lon1r = lon1 * d2r
  brngr = brng * d2r

  ; Calculate new latitude and longitude in radians
  lat2r = asin(sin(lat1r) * cos(d / R) + cos(lat1r) * sin(d / R) * cos(brngr))
  lon2r = lon1r + atan2( sin(brngr) * sin(d / R) * cos(lat1r), cos(d / R) - sin(lat1r) * sin(lat2r))

  ; Convert new latitude and longitude back to degrees
  lat2 = lat2r * r2d
  lon2 = lon2r * r2d

  ; Print the new latitude and longitude
  if (extrap_latlon) then
    print("Xtrap lat/lon:  " + lat2 +"/"+lon2)
  end if

  return_var = new(2,typeof(lat2))
  return_var(0) = lat2
  return_var(1) = lon2

  return(return_var)

end