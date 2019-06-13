;========================================================
; =============== COMMAND LINE ARGS
; YYYYMMMDDHH (string) --- time of CFSR analysis
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/storage/home/cmz5202/ncl/functions/namelist_extract.ncl"
;========================================================
begin

basehead = namelist_extract( nlfile, "basehead", "string")
fpname = namelist_extract( nlfile, "fpname", "string")
outdir = namelist_extract( nlfile, "outdir", "string")
wgt_file = namelist_extract( nlfile, "wgt_file", "string")

;f2name=fdir+"/"+"CHEY.VR28.NATL.EXT.CAM5.4CLM5.0.dtime900.cam.h2.1989-09-01-00000.nc"
;wgt_file="/global/homes/c/czarzyck/scratch/maps/hyperion/map_ne0np4natlanticext.ne30x4_to_0.25x0.25_GLOB.nc"
;outdir="/global/homes/c/czarzyck/scratch/et-hyperion/"
;basename=systemfunc("basename "+f2name)
;basename = str_sub_str(basename,"h2","h7")

;fpname="/global/homes/c/czarzyck/scratch/unigridFiles/ne0np4natlanticext.ne30x4/topo/topo_ne0np4natlanticext.ne30x4_smooth.nc"

f2=addfile(f2name,"r")
if (fpname .eq. "" .or. fpname .eq. "0") then
  fp="ZERO"
else
  fp=addfile(fpname,"r")
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
print(outdatestr)

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
tmp2d=f2->PSL(STTI:ENTI,:)
PS=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

tmp2d=f2->PSL(STTI:ENTI,:)
PSL=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

print("loading hybrid coefficients")
hyam=f2->hyam
hybm=f2->hybm
p0 = 1000.
;PRESSURE = pres_hybrid_ccm(PS,p0,hyam,hybm)

print("loading Z3")
tmp3d=f2->Z3(STTI:ENTI,:,:)
Z3=ESMF_regrid_with_weights(tmp3d,wgt_file,Opt)

;tmp2d=f2->T(STTI:ENTI,nlev-1,:)
;TBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

tmp2d=f2->UBOT(STTI:ENTI,:)
UBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

tmp2d=f2->VBOT(STTI:ENTI,:)
VBOT=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

tmp2d=f2->U850(STTI:ENTI,:)
U850=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

tmp2d=f2->V850(STTI:ENTI,:)
V850=ESMF_regrid_with_weights(tmp2d,wgt_file,Opt)

if (fp .eq. "ZERO") then
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
ilev = (/300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900./)
ilev!0 = "lev"
ilev&lev = ilev
ilev@units = "hPa"

intyp = 2      ; 1=linear, 2=log, 3=log-log
kxtrp = True   ; True=extrapolate
varflg = -1                          ; geo pot hgt is variable [tbot is used]


;Zi = vintp2p_ecmwf(Z3, PRESSURE/100., ilev ,PS, intyp, 1, \
;           kxtrp,varflg,TBOT,PHIS)
Zi = vinth2p(Z3,hyam,hybm,ilev,PS,intyp,p0/100.,1,kxtrp)


ZDIFF = Zi(:,{300.},:,:)
ZDIFF = Zi(:,{300.},:,:) - Zi(:,{500.},:,:)

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

end do

end

