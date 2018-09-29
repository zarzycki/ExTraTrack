undef("namelist_extract")
function namelist_extract(pfName[1]:string, keynam[1]:string, return_type[1]:string)
;--------------------------------------------------
;  Read a 'simple'  fortran namelist: key=value(s)
;  BEWARE: 'namelist_extract' has only been tested with one namelist. 
;  SUPPORT: None
;--------------------------------------------------
; Nomenclature:
;    pfName - path for file: dir+filname
;    keynam - key name 
;    return_type - "integer", "logical", "float", "string"
;
local namlst, nline, ngrp, n, chr, nchr, line, grpnam, tab, dq, delim \
     ,KEY, nKEY, iKEY, nfld
begin
   namlst  = asciiread(pfName, -1, "string") ; p(ath)/(f)ileName
   nline   = dimsizes(namlst)

;---Count the number of groups
   ngrp   = 0                        ; group counter
   do n=0,nline-1
      line = str_left_strip(namlst(n)) 
      chr := tochar(line)
      if (tostring(chr(0)).eq."&") then
          ngrp   = ngrp+1
          grpnam = tostring(chr(1:)) 
;;print("ngrp="+ngrp+"  group_name="+grpnam)
      end if
   end do

;---only one group allowed
   if (ngrp.gt.1) then
       print("namelist_extract:  ngrp="+ngrp)
       print("namelist_extract:  only one group allowed")
       exit
   end if

;--Look for line with 'keynam' to be parsed
   tab   = str_get_tab()
   dq    = str_get_dq()
   delim = "= ,"+tab            ; add 'tab'

   KEY   = str_left_strip(str_get_field(namlst, 1, delim))    ; extract leftmost field 
   nKEY  = dimsizes(KEY)
   iKEY  = ind(keynam.eq.KEY)                 ; index of line with 'keynam'

   if (ismissing(iKEY)) then
       print("namelist_extract:  key name="+keynam+" not found")
       print("namelist_extract:  KEY="+KEY)
       exit
   end if

;---Number of fields for 'namlst(iKEY)' line.  This includes the keyword.
   nfld  = str_fields_count(namlst(iKEY), delim)   

   chr   := tochar(namlst(iKEY))
   nchr   = dimsizes(chr)
   if (chr(nchr-1).eq.",") then
       nfld = nfld-1
   end if

   if (return_type.eq."string") then
       value=str_get_field(namlst(iKEY), 2, "=,"+tab+dq)
       print("Reading: "+keynam+" from "+pfName+" as "+return_type+":          "+value)
       return(value)
   end if

   ; if not string

   value = new(nfld, return_type)

   if (return_type.eq."integer") then
       do n=1,nfld
          value(n-1) = toint(str_get_field(namlst(iKEY), n+1, delim)) 
       end do
   end if

   if (return_type.eq."float") then
       do n=1,nfld
          value(n-1) = tofloat(str_get_field(namlst(iKEY), n+1, delim)) 
       end do
   end if

   if (return_type.eq."logical") then
       value = False                ; default value; NCL logical
       do n=1,nfld
          val = str_get_field(namlst(iKEY), n+1, delim)
          if (val .eq. "True" .or. val.eq.".true." .or. val.eq.".TRUE." .or. \
              val.eq.".T." .or. val.eq.".t." .or. val.eq."T") then
              value(n-1) = True     ; NCL logical
          end if   
       end do
   end if

   print("Reading: "+keynam+" from "+pfName+" as "+return_type+":          "+value)

   return(value) 
end


;;==================================================
;;==================================================
;;==================================================
;;&time_control
;; run_days                            = 5,
;; run_hours                           = 0,
;; run_minutes                         = 0,
;; run_seconds                         = 0,
;; start_year                          = 2015, 2001, 2001,
;; start_month                         = 06,   06,   06,
;; start_day                           = 19,   11,   11,
;; start_hour                          = 00,   12,   12,
;; start_minute                        = 00,   00,   00,
;; start_second                        = 00,   00,   00,
;; end_year                            = 2015, 2001, 2001,
;; end_month                           = 06,   06,   06,
;; end_day                             = 24,   12,   12,
;; end_hour                            = 00,   12,   12,
;; end_minute                          = 00,   00,   00,
;; end_second                          = 00,   00,   00,
;; interval_seconds                    = 600
;; input_from_file                     = .true.,.true.,.true.,
;; history_interval                    = 10,   60,   60,
;; frames_per_outfile                  = 1,    1,    1,
;; restart                             = .false.,
;; restart_interval                    = 5000,
;; io_form_history                     = 2
;; io_form_restart                     = 2
;; io_form_input                       = 2
;; io_form_boundary                    = 2
;; auxinput1_inname                    = "../../../WPS/wpsout_cus_4km/met_em.d<domain>.<date>"
;; history_outname                     = "wrfout/wrfout_d<domain>_<date>"
;; debug_level                         = 0
