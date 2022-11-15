;========================================================
; This file reads in a traj file produced by the TC tracker
; and produces a lat-lon plot of trajectories color-coded by intensity
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/glade/u/home/zarzycki/ncl/functions/surfaceWindCorrFactor.ncl"
load "/glade/u/home/zarzycki/ncl/functions/getTrajectoriesGeneric.ncl"
load "/glade/u/home/zarzycki/ncl/functions/mask_tc.ncl"

;========================================================
begin

thefile="./et-tracker/text_files/traj_et_cfsr_avg"
ncols = -1                    ; assign number of columns
maxLines = -1
headerDelimStr = "start"
delimStr="tab"
isHeader=False
traj_data = getTrajectoriesGeneric(thefile,ncols,maxLines,headerDelimStr,delimStr,False,isHeader)

lon = traj_data(:,0,:)
lat = traj_data(:,1,:)
wind = traj_data(:,3,:)
pres = traj_data(:,2,:)


B = traj_data(:,6,:)
Vlt = traj_data(:,7,:)
Vut = traj_data(:,8,:)

yyyy=  toint(traj_data(:,9,:))
mm= toint(traj_data(:,10,:))
dd= toint(traj_data(:,11,:))
hh= toint(traj_data(:,12,:))

ndims = dimsizes(lon)
nstorms = ndims(0)
ntimes = ndims(1)

print(nstorms+" "+ntimes)
doyarray=new(367,"integer")
doyarray=0
do ii = 0,nstorms-1
do jj = 0,ntimes-1
  if (.not. ismissing(B(ii,jj))) then
    if (Vlt(ii,jj) .gt. 0. .and. Vut(ii,jj) .gt. 0.) then
      thisDoy=day_of_year(yyyy(ii,jj),mm(ii,jj),dd(ii,jj))
      doyarray(thisDoy)=doyarray(thisDoy)+1
    end if
  end if
end do
end do

 wks   = gsn_open_wks ("x11","xy")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Basic XY plot"          ; add title

 plot  = gsn_csm_xy (wks,ispan(0,366,1),doyarray,res) ; create plot



exit

pdf2_nc  = pdfxy(ndtooned(Vlt),ndtooned(B),0,0,False)
pdf2_cn  = pdfxy(ndtooned(Vlt),ndtooned(Vut),0,0,False)


 wks      = gsn_open_wks("x11","pdf")           ; send graphics to PNG file

  res      = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@cnInfoLabelOn = False

  res@cnFillOn            = True    ; Turn on color
  res@cnFillPalette       = "amwg"  ; set color map
  res@cnLinesOn           = False
  res@cnLineLabelsOn      = False

  plot = gsn_csm_contour (wks,pdf2_nc, res)

  draw(plot)
  frame(wks)

  plot = gsn_csm_contour (wks,pdf2_cn, res)
  draw(plot)
  frame(wks)


end
