@all_tables.pro
pro impurity_table,qimp=qimp,plot=plot,ps=ps
  nmax=12                 ;; number of quantum states
  nmax_red=6              ;; reduced number of quantum states whereby the 
                          ;; states upto nmax are taken into account 
                          ;; as loss mechanism!
  maxeb=400.d0  ;; maximum energy in table [kev]
  maxti=20.d0   ;; maximum ion temperature in table [kev]         
  neb=201       ;; table size
  nti=51        ;; table size


  ;; arrays of eb and ti
  ebarr=maxeb*dindgen(neb)/(neb-1.)
  tiarr=maxti*dindgen(nti)/(nti-1.)
  deb=ebarr[2]-ebarr[1]
  dti=tiarr[2]-tiarr[1]
  print,'d-energy:', deb,'[kev]'
  print,'d-ti:'    , dti,'[kev]' 
  
  if not(keyword_set(qimp)) and not(keyword_set(inputs))then begin
     print, 'which imurity charge?'
     qimp=5
     read,qimp
  endif
  
  case qimp of 
     7: impchar='nitrogen'
     6: impchar='carbon'
     5: impchar='boron'
     else: stop
  endcase
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


     ;; CHARGE EXCHANGE 
     xtit='Energy [keV/amu]'
     ytit='Cross-section [cm!u2!n]'
     xran=[1.e-2,1.e3]
     yran=[1.e-20,1.e-12]    
 
     file='PLOTS/'+impchar+'_charge_exchange.eps'
     if keyword_set(ps) then device,filename=file else window,0
     plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
          ,xtit=xtit,ytit=ytit,/xstyle
     for i=0,13 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],col=0,linestyle=1
     for i=0,7 do oplot,[xran[0]*10.^i,xran[0]*10.^i],yran,color=0,linestyle=1
     for n=0, nmax-1 do begin ;; lower state
        oplot,ebarr,icx_adas_janev(ebarr,n+1,qimp),color=colorarr[n],thick=2
        xyouts,0.24,0.79-0.05*n,'n='+string(n+1,f='(1i2)') $
               ,color=colorarr[n],/norm
     endfor    
     xyouts,.22,0.86,'Charge-exchange reactions with '+impchar+' from: ',color=0,/norm 
  

 ;; IMPACT IONIZATOIN
     if keyword_set(ps) then device,file='PLOTS/'+impchar+'_impact_ionization.eps' else window,3
     plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
          ,xtit=xtit,ytit=ytit,/xstyle
     for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],col=0,linestyle=1
     for i=0,7 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
     for n=0,nmax-1 do begin ;; lower state
        oplot,ebarr,iiionization(ebarr,n+1,qimp),color=colorarr[n]
        xyouts,0.24,0.79- n*0.05,'n='+string(n+1,f='(1i3)'),/norm,color=colorarr[n]
     endfor
     xyouts,.22,0.86,impchar+' impact-ionization from: ',color=0,/norm 
  

     file='PLOTS/'+impchar+'_impact_excitation'
     if not keyword_set(ps) then window,1,xsize=1200,ysize=600
     for n=0,nmax-2 do begin         ;; lower state
        if keyword_set(ps) then device,filename=file+strtrim(string(n+1,f='(i2)'),2)+'.eps'
        plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
             ,xtit=xtit,ytit=ytit
        for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],color=0,linesty=1
        for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
   xyouts,.22,0.86,impchar+' impact-excitation from n='+string(n+1,f='(1i2)')+' to:',color=0,/norm
        for m=n+1,nmax-1 do begin    ;; upper state
           sigma=iiexcitation(ebarr,n+1,m+1,qimp)
           oplot,ebarr,sigma,color=colorarr[m]
           if total(sigma[where(finite(sigma) ne 0)]) gt 0 then begin
            xyouts,0.24,0.83- m*0.05, 'm='+string(m+1,f='(1i2)'),/norm, color=colorarr[m] 
           endif
        endfor
     endfor  


     if keyword_set(ps) then begin
        device, /close
        set_plot,'X'
        device,decomposed=0
     endif
     stop
  endif


  if qimp eq 5 then begin 
     aimp=10.8
     file='qbtable.bin'
     file_full='qbtable_full.bin'
  endif
  if qimp eq 6 then begin 
     aimp=12.0
     file='qctable.bin'
     file_full='qctable_full.bin'
  endif
  if qimp eq 7 then begin 
     aimp=14.0
     file='qntable.bin'
     file_full='qntable_full.bin'
  endif
  ab = 2.
  qi=fltarr(nmax+1,nmax,neb,nti) 
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
              qi[m,n,ie,iti]=beam_therm_rate(ti,eb,aimp,ab,de,'iiexcitation', $
                                             param=[n+1,m+1,qimp])
              qi[n,m,ie,iti]=beam_therm_rate(ti,eb,aimp,ab,de,'iiexcitation', $
                                             param=[n+1,m+1,qimp],/deexcit)
           endfor
        endfor 
        ;; Impact ionization ;; Charge exchange
        for n=0,nmax-1 do begin 
           de = 13.6d-3/(n+1)^2
           qi[nmax,n,ie,iti]= $
             beam_therm_rate(ti,eb,aimp,ab,de,'iiionization',param=[n+1,qimp]) $
           +beam_therm_rate(ti,eb,aimp,ab,0.0,'icx_adas_janev',param=[n+1,qimp])
        endfor
     endfor
  endfor
  ;; ---------------------------------------
  ;; ----- STORE DATE INTO BINARY FILES ----
  ;; ---------------------------------------
  openw, lun, file_full, /get_lun
  writeu,lun, long(nti)
  writeu,lun, double(dti)
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax) 
  writeu,lun, double(qi)
  close,lun
  free_lun, lun
  ;; ------- STORE DATE WITH REDUCED N-LEVELS ----
  reduce_table,qi,neb,nti,nmax_red,qi_red
  openw, lun, file, /get_lun
  writeu,lun, long(nti)
  writeu,lun, double(dti)
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax_red) 
  writeu,lun, double(qi_red)
  close,lun
  free_lun, lun
  print, 'impurity table written to:', file
end
