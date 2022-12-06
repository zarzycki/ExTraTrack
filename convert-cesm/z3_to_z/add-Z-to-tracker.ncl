begin

fileZ3=addfile(filename,"w")

; create p levels to interp to
pnew = (/300.,350.,400.,450.,500.,550.,600., \
                         650.,700.,750.,800.,850.,900./)

; load data
Z3=fileZ3->Z3
PS=fileZ3->PS
hyam=fileZ3->hyam
hybm=fileZ3->hybm
p0=1000.0

; interp Z3 to constant pnew levs
Z = vinth2p(Z3,hyam,hybm,pnew,PS,2,p0,1,True)

; write Z to original file
fileZ3->Z=Z

; NOTE: vinth2p seems to store coordinate variable as lev_p, and NCL will write this
; as a variable to NetCDF implicitly when adding Z, so let's comment this out for now
; so we don't get any warnings...
;lev_p=Z&lev_p
;fileZ3->lev_p=lev_p

end
