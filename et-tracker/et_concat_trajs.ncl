;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "../functions/namelist_extract.ncl"
;========================================================
; concat 
;========================================================
begin
print(" ")

etfile = namelist_extract( nlfile, "etfile", "string")
type = namelist_extract( nlfile, "type", "string")

purgetmpfiles=False

textfiledir="./text_files/"
if ( toint(systemfunc("test -d "+textfiledir+"; echo $?")) .ne. 0 ) then
  system("mkdir -p "+textfiledir+"/")
end if

system("cat ./text_files/tmp_"+type+"_*.txt > "+etfile)
if (purgetmpfiles) then
  system("tmp_"+type+"_*.txt")
end if
	
end
