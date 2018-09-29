undef("calc_B")
function calc_B(
  Z_thick : numeric,  
  qq_ddeg[1] : numeric,
  critRad[1] : numeric,
  dist[1] : numeric,
  ang[1] : numeric,
  minPlat[1] : numeric,
  minPlon[1] : numeric,
  opt[1] : numeric
  )
local pi, deg2rad, qq_lat, qq_lon, qq_nlat, qq_nlon, qq_ang_all, qq_wgt, Zl, Zr, d, dist, ang, \
  ndims, nlat, nlon, yyyy, mm, dd, hh, units, time, splitDate, dtime_map, outName, Zltmp, Zrtmp
begin

  pi = 4.*atan(1.)
  deg2rad = pi/180.

    ; assign left and right Zs to entire Z array
    Zl = Z_thick({(minPlat-qq_ddeg):(minPlat+qq_ddeg)},{(minPlon-qq_ddeg):(minPlon+qq_ddeg)})
    Zr = Zl
    
    ; subselect lat/lon arrays to within 10deg of minPlon
    qq_lat=Zl&lat
    qq_lon=Zl&lon
    qq_nlat = dimsizes(qq_lat)
    qq_nlon = dimsizes(qq_lon)
    qq_ang_all = new((/qq_nlat,qq_nlon/),float) ; angle of all values from storm center
	  qq_wgt = cos(qq_lat*deg2rad)
    
    do jlat = 0,qq_nlat-1
      do jlon = 0,qq_nlon-1
      thisLat=qq_lat(jlat)
      thisLon=qq_lon(jlon)
    
      ; find bearing between minimum pressure and each lat/lon pair
      qq_ang_all(jlat,jlon) = bearing(minPlat,minPlon,qq_lat(jlat),qq_lon(jlon))

      ; find 500 km radius for Z calculations	  
      d = gc_latlon(minPlat,minPlon,qq_lat(jlat),qq_lon(jlon),2,4)  ; km
      if (d .gt. critRad) then
        Zl(jlat,jlon) = Zl@_FillValue  ; set points outside radius to missing
        Zr(jlat,jlon) = Zr@_FillValue  ; set points outside radius to missing
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

    if (1 .eq. 0) then
      print("... writing out B to tmp netCDF")
      dtime_map = (/4,2,2,2/)
      splitDate = str_split_by_length(tostring(opt@date),dtime_map)
      yyyy=toint(splitDate(0))
      mm=toint(splitDate(1))
      dd=toint(splitDate(2))
      hh=toint(splitDate(3))
      units  = "days since 1900-01-01 00:00:00"
      time   = cd_inv_calendar(yyyy,mm,dd,hh,0,0,units,0)

      ndims=dimsizes(Zl)
      nlat=ndims(0)
      nlon=ndims(1)
      Zltmp=new((/1,nlat,nlon/),typeof(Zl))
      Zrtmp=new((/1,nlat,nlon/),typeof(Zr))

      Zltmp!0="time"
      Zltmp!1="lat"
      Zltmp!2="lon"
      Zltmp&time=time
      Zltmp&lat=Zl&lat
      Zltmp&lon=Zl&lon
      copy_VarMeta(Zltmp,Zrtmp)

      Zltmp(0,:,:)=Zl
      Zrtmp(0,:,:)=Zr

      outName="SAMP_"+tostring(opt@date)+"_"+tostring(opt@id)+".nc"
      system("/bin/rm -f "+outName)
      ncdf = addfile(outName ,"c")
      fAtt               = True
      fAtt@Conventions   = "None"
      fAtt@creation_date = systemfunc ("date")
      fileattdef( ncdf, fAtt )            ; copy file attributes
      filedimdef(ncdf,"time",-1,True) 
      ncdf->Zl=Zltmp
      ncdf->Zr=Zrtmp

      delete([/yyyy,mm,dd,hh,units,time,splitDate,dtime_map/])
      delete([/ndims,nlat,nlon/])
      delete([/outName,Zltmp,Zrtmp/])
    end if
  
    ;========================================================
    ; Compute B value
    ;========================================================

    ; cyclone thermal symmetry parameter (Eq. 2, Hart 2003) using weighted average
    B = wgt_areaave(Zr,qq_wgt,1.0,0) - wgt_areaave(Zl,qq_wgt,1.0,0)
  
    ; if southern hemisphere, h needs to be negative per Evans and Hart.
    ;Evans, J.L. and R.E. Hart, 2003: Objective Indicators of the Life Cycle Evolution
    ;of Extratropical Transition for Atlantic Tropical Cyclones. Mon. Wea. Rev., 131, 909–925,
    if (minPlat .lt. 0)
      B = -B
    end if

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
  indices[3] : numeric,
  logp[*] : numeric,
  opt[1] : numeric
  )
local qq_Z, d, dZ, indtop, indmid, indbot, qq_lat, qq_lon, qq_nlat, qq_nlon, Vlt, Vut

begin

; assign indices to be used later
indtop = indices(0)      ; top index for VTU
indmid = indices(1)     ; bottom index for VTU, top index for B, VTL
indbot = indices(2)     ; bottom index for B, VTL

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

;print(logp(indmid:indbot))
;print(dZ(indmid:indbot))

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




undef("get_sizeRad_from_lat")
function get_sizeRad_from_lat(
  minPlat[1] : numeric
)
local radlatcuts,radsrccuts

begin
  ; search radius lat thresholds
  radlatcuts=(/20.0 ,-99.9, 45.0/)
  radsrccuts=(/300.0,500.0,700.0/)
  if ( minPlat .lt. radlatcuts(0) ) then
    rad_tmp=radsrccuts(0)
  else if ( minPlat .ge. minPlat ) then
    rad_tmp=radsrccuts(2)
  else
    rad_tmp=radsrccuts(1)
  end if
  end if
  return(rad_tmp)
end


