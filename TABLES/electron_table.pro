@all_tables.pro
pro electron_table,plot=plot,ps=ps
  nmax=12                      ;; number of quantum states
  nmax_red=6                   ;; reduced number of quantum states whereby the 
                               ;; states upto nmax are taken into account 
                               ;; as loss mechanism!

  maxeb=400.d0  ;; maximum energy in table [kev]
  maxte=20.d0   ;; maximum electron temperature in table [kev]         
  neb=201       ;; table size
  nte=201       ;; table size

  ;; arrays of eb and ti
  ebarr=maxeb*dindgen(neb)/(neb-1.)
  tearr=maxte*dindgen(nte)/(nte-1.)
  deb=ebarr[2]-ebarr[1]
  dte=tearr[2]-tearr[1]
  print,'d-energy:',deb,'[kev]'
  print,'d-te    :'    , dte,'[kev]' 

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
     
     xran=[1.e-2,1.e3]
     yran=[1.e-20,1.e-12]
     xtit='Energy [keV]'
     ytit='Cross-section [cm!u2!n]'
   
     file='PLOTS/electron_impact_excitation_n='
     if not keyword_set(ps) then  window,0,xsize=1200,ysize=600
     !P.multi=[0,0,0]
     for n=0,nmax-2 do begin         ;; lower state
        if keyword_set(ps) then device,filename=file+strtrim(string(n+1,f='(i2)'),2) +'.eps'
        plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
             ,xtit=xtit,ytit=ytit
        for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],color=0,linesty=1
        for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
        xyouts,.3,0.85,'Electron impact-excitation from n='+string(n+1,f='(1i2)')+' to:',color=0,/norm
        for m=n+1,nmax-1 do begin    ;; upper state
           oplot,tearr,eiexcitation_janev2004(tearr,n+1,m+1),color=colorarr[m] 
           xyouts,0.83,0.85- m*0.05, 'm='+string(m+1,f='(1i2)'),/norm, color=colorarr[m] 
        endfor
     endfor   
 

     !P.multi=0
     file='PLOTS/electron_impact_ionization.eps'
     if keyword_set(ps) then device,filename=file,xsize=8,ysize=6 else window,1
     plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
          ,xtit=xtit,ytit=ytit
     for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i],col=0,linestyle=1
     for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran,color=0,linestyle=1
   
     for n=0,nmax-1 do begin ;; lower state
        sigma=eiionization_janev2004(tearr,n+1)
        oplot,tearr,sigma,color=colorarr[n]
         xyouts,0.83,0.85- n*0.05,'n='+string(n+1,f='(1i2)') $ 
               ,/norm,color=colorarr[n]
      endfor
  xyouts,.3,0.85,'Electron impact-ionization from: ',color=0,/norm
     if keyword_set(ps) then begin
        device, /close
        set_plot,'X'
        device,decomposed=0
     endif
     stop
  endif
  
  ab = 2.
  ame = 9.109d-31/1.661d-27  
  qe=fltarr(nmax+1,nmax,neb,nte) 
  for ie=0L,neb-1 do begin
     print, ie, ' of ',neb-1
     eb=ebarr[ie]
     for ite=0L,nte-1 do begin  
        te=tearr[ite]
        ;; Excitation and Deexcitation
        for n=0,nmax-2 do begin
           en = 13.6d-3*(1.-1./(n+1)^2)
           for m=n+1,nmax-1 do begin
              em=13.6d-3*(1.-1./(m+1)^2)
              de = em-en
              qe[m,n,ie,ite] = beam_therm_rate(te,eb,ame,ab,de,'eiexcitation_janev2004',param=[n+1,m+1])
              qe[n,m,ie,ite] = beam_therm_rate(te,eb,ame,ab,de,'eiexcitation_janev2004',param=[n+1,m+1],/deexcit)
           endfor
        endfor    
        ;; Impact ionization ;; Charge exchange
        for n=0,nmax-1 do begin 
           de = 13.6d-3/(n+1)^2
           qe[nmax,n,ie,ite]= $
           beam_therm_rate(te,eb,ame,ab,de,'eiionization_janev2004',param=[n+1])
        endfor
     endfor
  endfor

  ;; ---------------------------------------
  ;; ----- STORE DATE INTO BINARY FILES ----
  ;; ---------------------------------------
  file='qetable_full.bin'
  openw, lun, file, /get_lun
  writeu,lun, long(nte)
  writeu,lun, double(dte)
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax) 
  writeu,lun, double(qe)
  close,lun
  free_lun, lun
  ;; ------- STORE DATE WITH REDUCED N-LEVELS ----
  reduce_table,qe,neb,nte,nmax_red,qe_red
  file='qetable.bin'
  openw, lun, file, /get_lun
  writeu,lun, long(nte)
  writeu,lun, double(dte)
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax_red) 
  writeu,lun, double(qe_red)
  close,lun
  free_lun, lun
  print,'electron table written!'
end
