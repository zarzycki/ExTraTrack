function buildTimeLookup(filelist)
; Pass in a filelist of CAM files, build a lookup table for filename and time indices
local ftmp, nfiles, nftimes, timeArr, listArr, ii
begin
  ; Here we build a "lookup" table to find what times are on what files...
  ; get list of h3 files
  ftmp = addfile(filelist(0),"r")
  nfiles = dimsizes(filelist)
  nftimes = dimsizes(ftmp->time)
  timeArr = new((/nfiles,nftimes/),typeof(ftmp->time))
  listArr = new((/nfiles/),"string")
  do ii = 0,nfiles-1
    print("loading time for file #"+ii+" of "+(nfiles-1)+"  "+filelist(ii))
    ftmp = addfile(filelist(ii),"r")
    nftimes=dimsizes(ftmp->time)   ; needed since files at end of series may have less than rest
    timeArr(ii,0:nftimes-1) = ftmp->time
    listArr(ii)=filelist(ii)
  end do
  timeArr!0="files"
  timeArr&files=listArr

  timeunits=ftmp->time@units
  opt=0
  opt@calendar=ftmp->time@calendar

  timeArr@units=timeunits
  timeArr@calendar=ftmp->time@calendar
  delete(ftmp)
  return(timeArr)
end



function getTimeLookup(filelist,stashedfile)
local f
begin

  f = addfile(stashedfile,"r")

  timeArr = f->timeArr
  timeArr&files=filelist

  delete(f)
  return(timeArr)

end



function getFileAndTimeFromLookup(timeArr,thisTime)
local timeArr, thisTime, time1d, dsizes_time, indices, fileix, trackindex, needed_file, timediff,eps
begin

  ; small number for checking time sync issues
  eps=1.0e-6

  ; look up time ix
  time1d      = ndtooned(timeArr)
  dsizes_time = dimsizes(timeArr)

  ; use this to get *exact match* for thisTime to some time in the time1d array
  ;timediff=0.0
  ;indices  = ind_resolve(ind(time1d.eq.thisTime),dsizes_time)

  ; here, find the index that is *nearest* to the thisTime
  ; which allows for slightly out of sync dates (e.g., 1 hour)
  timediff=abs(time1d-thisTime)
  indices  = ind_resolve(minind(timediff),dsizes_time)
  if (min(timediff) .gt. eps) then
    print("WARNING: diff b/w best matched time in nc files & time from traj file = "+min(timediff))
    print("WARNING: thisTime is: "+cd_calendar(thisTime,3))
    ;print(time1d+" "+timediff)
    print("---> this may mean your times are out of sync!")
    return_error = toint(min(timediff))
  else
    return_error = 0
  end if

  if (any(ismissing(indices))) then
    trackindex=-1
    needed_file=""
  else
    fileix     = indices(0,0)
    trackindex = indices(0,1)
    needed_file=timeArr&files(fileix)
  end if
  return( [/ needed_file,trackindex,return_error /] )
end

function getCoords(filelist,latBounds)
local ftmp,lat,lon,timeunits,opt,latBounds,minLatOut,maxLatOut,minLonOut,maxLonOut
begin

  ; unpack latBounds var
  minLatOut=latBounds(0)
  maxLatOut=latBounds(1)
  minLonOut=latBounds(2)
  maxLonOut=latBounds(3)

  ; we only need first file since we can get them from all files
  ftmp = addfile(filelist(0),"r")
  lat=ftmp->lat({minLatOut:maxLatOut})
  lon=ftmp->lon({minLonOut:maxLonOut})
  timeunits=ftmp->time@units
  opt=0
  opt@calendar=ftmp->time@calendar
  delete(ftmp)
  return( [/ lat,lon,timeunits,opt /] )
end



function packLatLonBounds(MINLAT,MAXLAT,MINLON,MAXLON)
local MINLAT,MAXLAT,MINLON,MAXLON
begin
  tmp = new(4,typeof(MINLAT))
  tmp(0)=MINLAT
  tmp(1)=MAXLAT
  tmp(2)=MINLON
  tmp(3)=MAXLON
  return(tmp)
end



function setupArrays(nstorms,ntimes,lat,lon)
local nlat,nlon,lat2D,lon2D,varArray,nstorms,ntimes,lat,lon
begin
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  lat2D=conform_dims((/nlat,nlon/),lat,0)
  lon2D=conform_dims((/nlat,nlon/),lon,1)
  varArray=new((/nstorms,ntimes,nlat,nlon/),"float")

  varArray!0    = "storm"
  varArray!1    = "time"
  varArray!2    = "lat"
  varArray!3    = "lon"
  varArray&lat = lat
  varArray&lon = lon

  sttime=new(nstorms,"integer")
  entime=new(nstorms,"integer")
  durtime=new(nstorms,"integer")
  sttime!0    = "storm"
  entime!0    = "storm"
  durtime!0   = "storm"

  return( [/ lat2D,lon2D,varArray,sttime,entime,durtime /] )
end

function loopAndExtract(SNOWFALL,VARNAME,lat2D,lon2D,deltaDeg,nstorms,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,latBounds,offset)
local i,j,loaded_file,thisTime,thislat,thislon,retlist,needed_file,trackindex,maxLat,minLat,maxLon,minLon,gcdist,f1,VARNAME,lat2D,lon2D,deltaDeg,nstorms,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,minLatOut,maxLatOut,minLonOut,maxLonOut,offset
begin

  ; unpack latBounds var
  minLatOut=latBounds(0)
  maxLatOut=latBounds(1)
  minLonOut=latBounds(2)
  maxLonOut=latBounds(3)

  loaded_file=""   ; set loaded file to empty
  do j = 0,nstorms-1
    print(" Storm #"+(j+1)+" of "+nstorms)
    durtime(j)=0
    do i = 0,ntimes-1
      if (.not.ismissing(traj_lon(j,i))) then

        print(i+"    "+ traj_yyyy(j,i)+" "+traj_mm(j,i)+" "+traj_dd(j,i)+" "+traj_hh(j,i) )
        thisTime=cd_inv_calendar(traj_yyyy(j,i),traj_mm(j,i),traj_dd(j,i),traj_hh(j,i),0.0,0.0,timeunits,opt)
        if (i .eq. 0) then
          sttime(j)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))
        end if
        durtime(j) = durtime(j)+6
        entime(j)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))

        retlist=getFileAndTimeFromLookup(timeArr,thisTime)
        needed_file=retlist[0]
        trackindex=retlist[1]
        delete(retlist)

        if (trackindex .ge. 0) then
          if (needed_file .ne. loaded_file) then
            ; need to load file
            print("Currently loading... "+needed_file)
            f1 = addfile(needed_file,"r")
            loaded_file=needed_file
          else
            print("existing loaded file... "+loaded_file)
          end if

          thislat=traj_lat(j,i)
          thislon=traj_lon(j,i)

          maxLat=thislat+deltaDeg
          minLat=thislat-deltaDeg
          maxLon=thislon+deltaDeg
          minLon=thislon-deltaDeg

          gcdist=gc_latlon(thislat,thislon,lat2D,lon2D,2,2)

          print(trackindex+"")
          SNOWFALL(j,i,:,:)=tofloat((f1->$VARNAME$(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
          SNOWFALL(j,i,:,:)=where(gcdist.lt.deltaDeg,SNOWFALL(j,i,:,:),SNOWFALL@_FillValue)
          SNOWFALL(j,i,:,:)=SNOWFALL(j,i,:,:)*offset
        end if
      end if
    end do
  end do

  return(SNOWFALL)

end


function loopAndExtractSingle(SNOWFALL,VARNAME,lat2D,lon2D,deltaDeg,stormix,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,latBounds,offset)
local i,j,loaded_file,thisTime,thislat,thislon,retlist,needed_file,trackindex,maxLat,minLat,maxLon,minLon,gcdist,f1,VARNAME,lat2D,lon2D,deltaDeg,nstorms,ntimes,sttime,entime,durtime,traj_lon,traj_lat,traj_yyyy,traj_mm,traj_dd,traj_hh,timeArr,timeunits,opt,minLatOut,maxLatOut,minLonOut,maxLonOut,offset
begin

  nstorms=1

  ; unpack latBounds var
  minLatOut=latBounds(0)
  maxLatOut=latBounds(1)
  minLonOut=latBounds(2)
  maxLonOut=latBounds(3)

  loaded_file=""   ; set loaded file to empty
  j=stormix

    print(" Storm #"+(j+1)+" of "+nstorms)
    durtime(0)=0
    do i = 0,ntimes-1
      if (.not.ismissing(traj_lon(j,i))) then

        print(i+"    "+ traj_yyyy(j,i)+" "+traj_mm(j,i)+" "+traj_dd(j,i)+" "+traj_hh(j,i) )
        thisTime=cd_inv_calendar(traj_yyyy(j,i),traj_mm(j,i),traj_dd(j,i),traj_hh(j,i),0.0,0.0,timeunits,opt)
        if (i .eq. 0) then
          sttime(0)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))
        end if
        durtime(0) = durtime(0)+6
        entime(0)=toint(traj_yyyy(j,i))*1000000+toint(traj_mm(j,i))*10000+toint(traj_dd(j,i))*100+toint(traj_hh(j,i))

        retlist=getFileAndTimeFromLookup(timeArr,thisTime)
        needed_file=retlist[0]
        trackindex=retlist[1]
        delete(retlist)

        if (trackindex .ge. 0) then
          if (needed_file .ne. loaded_file) then
            ; need to load file
            print("Currently loading... "+needed_file)
            f1 = addfile(needed_file,"r")
            loaded_file=needed_file
          else
            print("existing loaded file... "+loaded_file)
          end if

          thislat=traj_lat(j,i)
          thislon=traj_lon(j,i)

          maxLat=thislat+deltaDeg
          minLat=thislat-deltaDeg
          maxLon=thislon+deltaDeg
          minLon=thislon-deltaDeg

          gcdist=gc_latlon(thislat,thislon,lat2D,lon2D,2,2)

          print(trackindex+"")
          SNOWFALL(0,i,:,:)=tofloat((f1->$VARNAME$(trackindex,{minLatOut:maxLatOut},{minLonOut:maxLonOut})))
          SNOWFALL(0,i,:,:)=where(gcdist.lt.deltaDeg,SNOWFALL(0,i,:,:),SNOWFALL@_FillValue)
          SNOWFALL(0,i,:,:)=SNOWFALL(0,i,:,:)*offset
        end if
      end if
    end do

  return(SNOWFALL)

end



