undef("et_phase_plot")
procedure et_phase_plot(
  type[1]    : string,   ; type of data
  gifsdir[1] : string,   ; base directory to store images
  itc[1]     : integer,  ; storm number index
  date       : numeric,  ; date for plot title
  B          : numeric,  ; B parameter
  Vlt        : numeric,  ; lower thermal wind
  Vut        : numeric,  ; upper thermal wind
  debug      : logical   ; if true, optional debugging info
  )
local plot,wks,res,resP,txres,nom,dum0,dum1
; plots the sea level pressure and cyclone center
;
; D.R. Thatcher 23/07/2014
begin

  plot = new(2,graphic)
  plot2 = new(2,graphic)
  plot3 = new(2,graphic)

  system("mkdir -p "+gifsdir+"/phase_plots/")
  outfilename=gifsdir+"/phase_plots/phase_plot_"+type+"_"+sprinti("%0.3i",itc)
  if (debug) then
    print("Writing: "+outfilename)
  end if
  wks = gsn_open_wks("pdf",outfilename)

  res           = True
  res@gsnDraw   = False
  res@gsnFrame  = False

  res@tiMainString = " "
  res@gsnLeftString = "Start Date: "+date(0)
  res@gsnRightString = "End Date: "+date(dimsizes(date)-1)
  res@tmXBLabelFontHeightF   = 0.02
  res@tmYLLabelFontHeightF   = 0.02
  res@tiXAxisFontHeightF     = 0.025
  res@tiYAxisFontHeightF     = 0.025
  res@gsnXRefLine            = 0
  res@gsnYRefLine            = 10
  res@trXMinF                = -500
  res@trXMaxF                = 500

; First draw lines for both plots

  res@xyMarkLineMode         ="Lines"
  res@xyMarkerColor          = "Black"

  res@trYMinF                = -25
  res@trYMaxF                = 100
  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
  res@tiYAxisString          = "B [900-600 hPa Storm-Relative Thickness Symmetry]"

  plot(0)=gsn_csm_xy(wks,Vlt(::2),B(::2),res)

  res@trYMinF                = -500
  res@trYMaxF                = 500
  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
  res@tiYAxisString          = "-V~B~T~N~~S~U~N~ [600-300 hPa Thermal Wind]"

  plot(1)=gsn_csm_xy(wks,Vlt(::2),Vut(::2),res)

  delete(res@gsnLeftString)
  delete(res@gsnRightString)
  delete(res@gsnXRefLine)
  delete(res@gsnYRefLine)

; Overlay markers

  res@xyMarkLineMode         ="Markers"
  res@xyMarker               = 16
  res@xyMarkerSizeF          = 0.025
  res@xyMarkerColor          = "Grey50"

  res@trYMinF                = -25
  res@trYMaxF                = 100
  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
  res@tiYAxisString          = "B [900-600 hPa Storm-Relative Thickness Symmetry]"

  plot2(0)=gsn_csm_xy(wks,Vlt(::2),B(::2),res)
  overlay(plot(0),plot2(0))

  res@trYMinF                = -500
  res@trYMaxF                = 500
  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
  res@tiYAxisString          = "-V~B~T~N~~S~U~N~ [600-300 hPa Thermal Wind]"

  plot2(1)=gsn_csm_xy(wks,Vlt(::2),Vut(::2),res)
  overlay(plot(1),plot2(1))

; Overlay Marker Outline

;  res@xyMarkLineMode         ="Markers"
;  res@xyMarker               = 4
;  res@xyMarkerSizeF          = 0.025
;  res@xyMarkerColor          = "Black"

;  res@trYMinF                = -25
;  res@trYMaxF                = 100
;  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
;  res@tiYAxisString          = "B [900-600 hPa Storm-Relative Thickness Symmetry]"

;  plot3(0)=gsn_csm_xy(wks,Vlt(::2),B(::2),res)
;  overlay(plot(0),plot3(0))

;  res@trYMinF                = -500
;  res@trYMaxF                = 500
;  res@tiXAxisString          = "-V~B~T~N~~S~L~N~ [900-600 hPa Thermal Wind]"
;  res@tiYAxisString          = "-V~B~T~N~~S~U~N~ [600-300 hPa Thermal Wind]"

;  plot3(1)=gsn_csm_xy(wks,Vlt(::2),Vut(::2),res)
;  overlay(plot(1),plot3(1))

; add text

  txres               = True
  txres@txFontHeightF = 0.025
  txres@txFont        = "helvetica-bold"
  txres@txFontColor   = (/"Red3"/)

  nom=ispan(1,dimsizes(Vlt(::4)),1)
  dum0=gsn_add_text(wks,plot(0),sprintf("%6.4g", nom), Vlt(::4),B(::4),txres)
  dum1=gsn_add_text(wks,plot(1),sprintf("%6.4g", nom), Vlt(::4),Vut(::4),txres)

; Panel plots

  resP = True
  resP@txString = "Storm Number: "+itc
  gsn_panel(wks,plot,(/2,1/),resP)

  delete([/wks,plot,plot2,plot3,res,resP,txres,nom,dum0,dum1/])

end
