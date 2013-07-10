@all_tables.pro
pro proton_table,low_res=low_res,plot=plot,ps=ps
  print,'Calcuation of the proton table!'

  nmax=12                    ; number of quantum states
  maxeb=1000.d0              ;[kev]
  maxti=20.d0                ;[kev]

  neb=10001
  nti=201
  if keyword_set(low_res) then begin
     neb=101
     nti=21
  endif
  
  ;; arrays of eb and ti
  ebarr=maxeb*dindgen(neb)/(neb-1.)
  tiarr=maxti*dindgen(nti)/(nti-1.)
  deb=ebarr[2]-ebarr[1]
  dti=tiarr[2]-tiarr[1]
  print,'d-energy:', deb,'[kev]'
  print,'d-ti:    '    , dti,'[kev]' 
  ;; plot cross sections
  if keyword_set(plot) then begin
     loadct,39         ;;green=150;blue=50;yellow=200;red=254;black=0;white=255
     tvlct,0,150,0,150 ;; darker green!
     linthick=1.
     set_plot,'X' & device, decomposed=0
     !P.charthick=1. & !P.charsize=1.5 & !P.thick=1.0
     if keyword_set(ps) then begin  
        set_plot, 'ps'
        device, font_size=8, inches=0 , /encaps $
                ,xsize=7.9, ysize=6, /color, bits_per_pixel=8 
        device, /helvetica,font_index=3
        device, /symbol,font_index=4                                    
        device, /helvetica
        linthick=1.5
        !P.charsize=1. & !P.charthick=3. & !P.thick=1.5 & !P.font=0
     endif
     !P.background=255 & !P.color=0 &  !p.multi=0 
     colorarr=[0,20,40,60,80,100,120,140,160,180,200,220]

     file='PLOTS/proton_impact_ionization.eps'
     if keyword_set(ps) then device,filename=file else window,2
     xtit='Energy [keV/amu]'
     ytit='Cross-section [cm!u2!n]'
     xran=[1.e-2,1.e3]
     yran=[1.e-20,1.e-12]
     plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
          ,xtit=xtit,ytit=ytit
     for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],color=0,linesty=1
     for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
     for n=0,nmax-1 do begin ;; lower state
        sigma=omullane(ebarr,n+1)
        oplot,ebarr,sigma,color=colorarr[n]
        xyouts,0.24,0.79- n*0.05,'n='+string(n+1,f='(1i2)'),/norm,color=colorarr[n]
     endfor
     xyouts,.22,0.86,'Proton impact-ionization from: ',color=0,/norm 


     if not keyword_set(ps) then window,1,xsize=1200,ysize=600
     file='PLOTS/proton_impact_excitation'
     for n=0,nmax-2 do begin         ;; lower state
        if keyword_set(ps) then device,filename=file+strtrim(string(n+1,f='(i2)'),2)+'.eps'
        plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
             ,xtit=xtit,ytit=ytit
        for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],color=0,linesty=1
        for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
        xyouts,.22,0.86,'Proton impact-excitation from n='+string(n+1,f='(1i2)')+' to:',color=0,/norm
        for m=n+1,nmax-1 do begin    ;; upper state
           sigma=piexcitation_janev(ebarr,n+1,m+1)
           oplot,ebarr,sigma,color=colorarr[m]
           xyouts,0.24,0.83- m*0.05, 'm='+string(m+1,f='(1i2)'),/norm, color=colorarr[m] 
        endfor
     endfor  

     file='PLOTS/proton_charge_exchange.eps'
     if keyword_set(ps) then device,filename=file else window,0
     plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
          ,xtit=xtit,ytit=ytit
     for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],color=0,linesty=1
     for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
     for n=0,nmax-1 do begin
        sigma=sigma_cx(ebarr,n+1,-1)
        oplot,ebarr,sigma,color=colorarr[n]
        xyouts,0.24,0.79-0.05*n,'n='+string(n+1,f='(1i2)'),/norm,color=colorarr[n]
     endfor
     xyouts,.22,0.86,'Charge-exchange reactions with protons from: ',color=0,/norm 


     if keyword_set(ps) then begin
        device, /close
        set_plot,'X'
        device,decomposed=0
     endif
     stop
  endif
  

  ab = 2.
  ai = 2.
  qp=fltarr(nmax+1,nmax,neb,nti) 
  for ie=0,neb-1 do begin
     if ie mod 10 eq 0 then print, ie, ' of ',neb-1
     eb=ebarr[ie]
     for iti=0L,nti-1 do begin  
        ti=tiarr[iti]
        ;; Excitation and Deexcitation
        for n=0,nmax-2 do begin
           en = 13.6d-3*(1.-1./(n+1)^2)
           for m=n+1,nmax-1 do begin
              em=13.6d-3*(1.-1./(m+1)^2)
              de = em-en
              qp[m,n,ie,iti] = beam_therm_rate(ti,eb,ai,ab,de,'piexcitation_janev',param=[n+1,m+1])
              qp[n,m,ie,iti] = beam_therm_rate(ti,eb,ai,ab,de,'piexcitation_janev',param=[n+1,m+1],/deexcit)
           endfor
        endfor    
        ;; Impact ionization ;; Charge exchange
        for n=0,nmax-1 do begin 
           de = 13.6d-3/(n+1)^2
           qp[nmax,n,ie,iti]= $
              beam_therm_rate(ti,eb,ai,ab,de $
                              ,'omullane',param=[n+1]) $
              + beam_therm_rate(ti,eb,ai,ab,0.0 $
                                ,'sigma_cx',param=[n+1,-1])
        endfor
     endfor
  endfor
  file='qptable.bin'
  openw, lun, file, /get_lun
  writeu,lun, long(nti)
  writeu,lun, double(dti)
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax) 
  for n=0,nmax-1 do begin
     for m=0,nmax do begin
        for ie=0,neb-1 do begin   
           for iti=0,nti-1 do begin       
              writeu, lun, float(qp[m,n,ie,iti])
           endfor
        endfor
     endfor
  endfor
  close,lun
  free_lun, lun
  print, 'proton table written to:', file
end
