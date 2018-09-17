undef("calc_B")
function calc_B(
  Z_thick : numeric,  
  lat[*] : numeric,
  lon[*] : numeric,
  qq_ddeg[1] : numeric,
  critRad[1] : numeric,
  dist[1] : numeric,
  ang[1] : numeric,
  minPlat[1] : numeric,
  minPlon[1] : numeric
  )
local pi, deg2rad, qq_lat, qq_lon, qq_nlat, qq_nlon, qq_ang_all, qq_wgt, Zl, Zr, d, dist, ang

begin

  pi = 4.*atan(1.)
  deg2rad = pi/180.
    
    ; subselect lat/lon arrays to within 10deg of minPlon
    qq_lat = lat({(minPlat-qq_ddeg):(minPlat+qq_ddeg)})
    qq_lon = lon({(minPlon-qq_ddeg):(minPlon+qq_ddeg)})
    qq_nlat = dimsizes(qq_lat)
    qq_nlon = dimsizes(qq_lon)
    qq_ang_all = new((/qq_nlat,qq_nlon/),float) ; angle of all values from storm center
	  qq_wgt = cos(qq_lat*deg2rad)

    ; assign left and right Zs to entire Z array
    Zl = Z_thick({(minPlat-qq_ddeg):(minPlat+qq_ddeg)},{(minPlon-qq_ddeg):(minPlon+qq_ddeg)})
    Zr = Zl
    
    ; find bearing between minimum pressure and each lat/lon pair
    do jlat = 0,qq_nlat-1
      do jlon = 0,qq_nlon-1
      thisLat=qq_lat(jlat)
      thisLon=qq_lon(jlon)
    
      ; only update angles if storm has moved at least 10 km
      if (dist .ge. 1 .or. any(ismissing(qq_ang_all(jlat,jlon)))) then
        qq_ang_all(jlat,jlon) = bearing(minPlat,minPlon,qq_lat(jlat),qq_lon(jlon))
      end if

      ; find 500 km radius for Z calculations	  
      d = gc_latlon(minPlat,minPlon,qq_lat(jlat),qq_lon(jlon),2,4)  ; km
      if (d .gt. critRad) then
        Zl(jlat,jlon) = Zl@_FillValue  ; set points outside radius to missing
        Zr(jlat,jlon) = Zr@_FillValue  ; set points outside radius to missing
        ;qq_Z(:,jlat,jlon) = qq_Z@_FillValue
      end if 
    
      ; set values along great circle line to missing (neither left nor right)
      if (qq_ang_all(jlat,jlon) .eq. ang) then
        Zl(jlat,jlon) = Zl@_FillValue
        Zr(jlat,jlon) = Zr@_FillValue
      
      ; for storm travel angles in quadrants 1 and 2 (NE and SE)
      else if (ang .ge. 0 .and. ang .lt. 180) then
        if (qq_ang_all(jlat,jlon) .gt. ang .and. qq_ang_all(jlat,jlon) .lt. ang+180) then
          Zl(jlat,jlon) = Zl@_FillValue
        else
          Zr(jlat,jlon) = Zr@_FillValue
        end if
      
      ; for storm travel angles in quadrants 3 and 4 (NW and SW)
      else if (ang .ge. 180 .and. ang .lt. 360) then
        if (qq_ang_all(jlat,jlon) .gt. ang-180 .and. qq_ang_all(jlat,jlon) .lt. ang) then
          Zr(jlat,jlon) = Zr@_FillValue
        else
          Zl(jlat,jlon) = Zl@_FillValue
        end if	  
      end if
      end if
      end if
    
      delete(d)
    
      end do
    end do
  
    delete([/jlat,jlon/])
  
    ;========================================================
    ; Compute B value
    ;========================================================

    ; cyclone thermal symmetry parameter (Eq. 2, Hart 2003) using weighted average
    B = wgt_areaave(Zr,qq_wgt,1.0,0) - wgt_areaave(Zl,qq_wgt,1.0,0)
  
    delete([/Zl,Zr/])

return(B)

end











undef("calc_VltVut")
function calc_VltVut(
  Z : numeric,  
  lat[*] : numeric,
  lon[*] : numeric,
  qq_ddeg[1] : numeric,
  critRad[1] : numeric,
  dist[1] : numeric,
  ang[1] : numeric,
  minPlat[1] : numeric,
  minPlon[1] : numeric,
  logp[*] : numeric
  )
local qq_Z, d, dZ, indtop, indmid, indbot, qq_lat, qq_lon, qq_nlat, qq_nlon, Vlt, Vut

begin

; assign indices to be used later
indtop = 0      ; top index for VTU
indmid = 6      ; bottom index for VTU, top index for B, VTL
indbot = 12     ; bottom index for B, VTL

qq_Z = Z(:,{(minPlat-qq_ddeg):(minPlat+qq_ddeg)},{(minPlon-qq_ddeg):(minPlon+qq_ddeg)})

qq_lat = lat({(minPlat-qq_ddeg):(minPlat+qq_ddeg)})
qq_lon = lon({(minPlon-qq_ddeg):(minPlon+qq_ddeg)})
qq_nlat = dimsizes(qq_lat)
qq_nlon = dimsizes(qq_lon)

do jlat = 0,qq_nlat-1
  do jlon = 0,qq_nlon-1
    ; find 500 km radius for Z calculations	  
    d = gc_latlon(minPlat,minPlon,qq_lat(jlat),qq_lon(jlon),2,4)  ; km
    if (d .gt. critRad) then
      qq_Z(:,jlat,jlon) = qq_Z@_FillValue
    end if 
  end do
end do

; find dZ = Zmax - Zmin at each pressure level
dZ = new((/dimsizes(logp)/),float)
do jlev = 0,dimsizes(logp)-1
  dZ(jlev) = max(qq_Z(jlev,:,:)) - min(qq_Z(jlev,:,:))
end do

; Compute thermal winds using Eqs. from Hart 2003 using linear regression
Vlt = regline(logp(indmid:indbot),dZ(indmid:indbot))  ; lower troposphere (Eq. 5)
Vut = regline(logp(indtop:indmid),dZ(indtop:indmid))  ; upper troposphere (Eq. 6)

returnvar=new(2,"float")
returnvar(0)=Vlt
returnvar(1)=Vut
return(returnvar)

end





undef("print_ET_toScreen")
procedure print_ET_toScreen(
  minP[1] : numeric,
  minPlat[1] : numeric,
  minPlon[1] : numeric,
  dist[1] : numeric,
  V[1] : numeric,
  ang[1] : numeric,
  d_ang[1] : numeric,
  B[1] : numeric,
  Vlt[1] : numeric,
  Vut[1] : numeric,
  TYPE[1] : string
)
begin
      print("   " + sprintf("%6.0f",minP) + \
            "   " + sprintf("%5.2f",minPlat) + \
            "   " + sprintf("%6.2f",minPlon) + \
            "   " + sprintf("%5.1f",dist) + \
            "   " + sprintf("%4.1f",V) + \
            "   " + sprintf("%5.1f",ang) + \
            "   " + sprintf("%6.1f",d_ang) + \
            "   " + sprintf("%5.1f",B) + \
            "   " + sprintf("%4.0f",Vlt) + \
            "   " + sprintf("%4.0f",Vut) + \
            "   " + "   "+TYPE)
end


