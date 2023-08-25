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

defaultnl="./namelist"
if (isdefined("nlfile")) then
  print("reading namelist from "+nlfile)
else if ((.not. isdefined("nlfile")) .and. (fileexists(defaultnl))) then
  print("using default namelist "+defaultnl)
  nlfile=defaultnl
else
  print("no namelist file defined, exiting...")
  print("feed in namelist file on the command line as in README")
  exit
end if
end if

etfile = namelist_extract( nlfile, "etfileori", "string")
type = namelist_extract( nlfile, "type", "string")

purgetmpfiles=True

textfiledir="./text_files/"
if ( toint(systemfunc("test -d "+textfiledir+"; echo $?")) .ne. 0 ) then
  system("mkdir -p "+textfiledir+"/")
end if

system("cat ./text_files/tmp_"+type+"_*.txt > "+etfile)
if (purgetmpfiles) then
  system("cd ./text_files/ ; rm -v tmp_"+type+"_*.txt")
end if

end
