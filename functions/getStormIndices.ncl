function getStormIndices(etfile)

; function takes in an ExTraTrack et_file and loads a 1-D integer array of storm indices
; which are the last column in the header row

local nstorms_tot, etdata_tmp, n1, nrows, dims

begin

nstorms_tot = stringtoint(systemfunc("grep 'start' "+etfile+" | wc -l"))
etdata_tmp = asciiread(etfile,-1,"string")  ; read in trajectory data as strings
dims = dimsizes(etdata_tmp)    ; get dimensions
nrows = dims(0)               ; get number of rows

s_num   = new(nstorms_tot,integer)

; loop over lines of trajectory data, convert string to float
nstorm = 0   ; assign storm number
do n1 = 0,nrows-1  ; loop over all data with header lines
  first = str_get_cols(etdata_tmp(n1),0,0)   ; obtain first character of row
  if (first .eq. "s") then              ; if not 's', this is not a header
    s_num(nstorm)  =  stringtoint(str_get_field(etdata_tmp(n1),7," "))
    nstorm = nstorm + 1
  end if
end do

delete([/nstorm,nstorms_tot,nrows,n1,etdata_tmp,dims/])

return(s_num)

end
