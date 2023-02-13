;========================================================
; =============== COMMAND LINE ARGS
; YYYYMMMDDHH (string) --- time of CFSR analysis
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/glade/u/home/zarzycki/sw/ExTraTrack/functions/namelist_extract.ncl"
;========================================================
begin

; 
; ncl create-files.ncl 'f2name="/glade/u/home/zarzycki/scratch/hyperion/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900/atm/CHEY.VR28.NATL.REF.CAM5.4CLM5.0.dtime900.cam.h2.2003-08-19-00000.nc"' 'nlfile="nl.hyperion"'     "

basehead = namelist_extract( nlfile, "basehead", "string")
fpname = namelist_extract( nlfile, "fpname", "string")
outdir = namelist_extract( nlfile, "outdir", "string")
wgt_file = namelist_extract( nlfile, "wgt_file", "string")
zleveloption = namelist_extract( nlfile, "zleveloption", "string")

;f2name=fdir+"/"+"CHEY.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.cam.h2.1989-09-01-00000.nc"
;wgt_file="/global/homes/c/czarzyck/scratch/maps/hyperion/map_ne0np4natlanticext.ne30x4_to_0.25x0.25_GLOB.nc"
;outdir="/global/homes/c/czarzyck/scratch/et-hyperion/"
;basename=systemfunc("basename "+f2name)
;basename = str_sub_str(basename,"h2","h7")
;fpname="/global/homes/c/czarzyck/scratch/unigridFiles/ne0np4natlanticext.ne30x4/topo/topo_ne0np4natlanticext.ne30x4_smooth.nc"

f2=addfile(f2name,"r")
if (fpname .eq. "" .or. fpname .eq. "0") then
  fpname="ZERO"
else
  fp=addfile(fpname,"r")
end if

; Special treatment for Hyperion because PSL instant lives on h3 files
special_treatment_hyperion=True
if (special_treatment_hyperion) then
  f3name = str_sub_str(f2name,".h2.",".h3.")
  f3=addfile(f3name,"r")
end if

nlev=dimsizes(f2->lev)
ntim=dimsizes(f2->time)

nloops = ntim/4  ; 6 hourly

do zz = 0,nloops-1

STTI=zz*4
ENTI=((zz+1)*4)-1

time=f2->time(STTI:ENTI)
utc_date=cd_calendar(time(0),0)
yyyy=toint(utc_date(:,0))
mm=toint(utc_date(:,1))
dd=toint(utc_date(:,2))
outdatestr=sprinti("%0.4i",yyyy)+"-"+sprinti("%0.2i",mm)+"-"+sprinti("%0.2i",dd)+"-00000"
basename=basehead+outdatestr+".nc"
outfile=outdir+"/"+basename
print(""+outdatestr)

;CORI.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.003.cam.h2.1984-03-31-00000.nc

;PSL=f2->PSL(0:ENTI,:)
;U850=f2->U850(0:ENTI,:)
;UBOT=f2->U(:,nlev-1,:)
;V850=f2->V850(0:ENTI,:)
;VBOT=f2->V(:,nlev-1,:)

Opt=True
;print("loading PRESSURE")
;tmp3d=f2->PRESSURE(STTI:ENTI,:,:)
;PRESSURE=ESMF_regrid_with_weights(tmp3d,wgt_file,Opt)

print("loading PS")
tmp2d=f2->PS(STTI:ENTI,:)
PS=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

print("loading PSL")
if (special_treatment_hyperion) then
  tmp2d=f3->PSL({time},:)
else
  tmp2d=f2->PSL(STTI:ENTI,:)
end if
PSL=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

print("loading hybrid coefficients")
hyam=f2->hyam
hybm=f2->hybm
p0 = 100000.
;PRESSURE = pres_hybrid_ccm(PS,p0,hyam,hybm)
intyp = 2      ; 1=linear, 2=log, 3=log-log
kxtrp = True   ; True=extrapolate
varflg = -1                          ; geo pot hgt is variable [tbot is used]

print("loading Z3")
tmp3d=f2->Z3(STTI:ENTI,:,:)
Z3=ESMF_regrid_with_weights(tmp3d,wgt_file,Opt)

;tmp2d=f2->T(STTI:ENTI,nlev-1,:)
;TBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

if (isfilevar(f2,"UBOT")) then
  tmp2d=f2->UBOT(STTI:ENTI,:)
else
  tmp2d=f2->U(STTI:ENTI,nlev-1,:)
end if
UBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)


if (isfilevar(f2,"VBOT")) then
  tmp2d=f2->VBOT(STTI:ENTI,:)
else
  tmp2d=f2->V(STTI:ENTI,nlev-1,:)
end if
VBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)


if (isfilevar(f2,"U850")) then
  tmp2d=f2->U850(STTI:ENTI,:)
else
  tmp3d=f2->U(STTI:ENTI,:,:)
  tmp2d=vinth2p(tmp3d,hyam,hybm,(/850./),PS,intyp,p0/100.,1,kxtrp)
end if
U850=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)


if (isfilevar(f2,"V850")) then
  tmp2d=f2->V850(STTI:ENTI,:)
else
  tmp3d=f2->V(STTI:ENTI,:,:)
  tmp2d=vinth2p(tmp3d,hyam,hybm,(/850./),PS,intyp,p0/100.,1,kxtrp)
end if
V850=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)


if (fpname .eq. "ZERO") then
  tmp1d=V850
  tmp1d=0.
  PHIS=tmp1d
else
  tmp1d=fp->PHIS(:)
  PHIS=ESMF_regrid_with_weights(tmp1d,wgt_file,Opt)
end if


delete(tmp1d)
delete(tmp2d)
delete(tmp3d)

;ilev = (/300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 775., 800., 825., 850., \
;  875., 900., 925., 950., 975., 1000/)
if (zleveloption .eq. "simple") then
  ilev = (/300., 500., 600., 900./)
else
  ilev = (/300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900./)
end if
ilev!0 = "lev"
ilev&lev = ilev
ilev@units = "hPa"

;Zi = vintp2p_ecmwf(Z3, PRESSURE/100., ilev ,PS, intyp, 1, \
;           kxtrp,varflg,TBOT,PHIS)
Zi = vinth2p(Z3,hyam,hybm,ilev,PS,intyp,p0/100.,1,kxtrp)

ZDIFF = Zi(:,{300.},:,:)
ZDIFF = Zi(:,{300.},:,:) - Zi(:,{500.},:,:)

system("mkdir -v -p "+outdir)
system("/bin/rm -f "+outfile)   ; remove any pre-existing file
ncdf = addfile(outfile,"c")  ; open output netCDF file

fAtt               = True
fAtt@source_file   =  f2name
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )
filedimdef(ncdf,"time",-1,True) 

ncdf->Z    = Zi
ncdf->Z300 = Zi(:,{300.},:,:)
ncdf->Z500 = Zi(:,{500.},:,:)
ncdf->DZ300500 = ZDIFF
ncdf->VBOT = VBOT
ncdf->UBOT = UBOT
ncdf->V850 = V850
ncdf->U850 = U850
ncdf->PSL  = PSL
ncdf->PS = PS

end do

end

