Function sigma1_astro,erel,m
  emin=5.
  emax=80.
  if m gt 4 then return,0.
  if erel lt emin or erel gt emax then return,0.
  x=(alog(erel/emin)-alog(emax/erel))/alog(emax/emin)  
  C=fltarr(8)
  C[1]=x
  C[2]=2.*x^2-1.
  C[3]=4.*x^3-3.*x
  C[4]=8.*(x^4-x^2)+1.
  C[5]=16.*x^5-20.*x^3+5.*x
  C[6]=32.*x^6-48.*x^4-18.*x^2-1.
  C[7]=64.*x^7-112.*x^5+56.*x^3-7.*x
  
  ;; FIT coefficients 
  ;; for n=1 to n=3:
  ;;1=s, 2=p, 3=d, 4=f
  D=fltarr(8,4)  
  E=fltarr(8,4)
  D[*,0]= [1.15,   0.963, -1.35, -0.242,  0.14,    0.00596,-0.0572, 0.0564  ]
  D[*,1]= [1.14,  -0.699, -1.05, -0.0922, 0.0511, -0.0818,-0.054,  -0.04    ]
  D[*,2]= [-1.65, -2.38,  -0.763, 0.137, -0.00275,-0.0614, 0.0162, -0.0629  ]
  ;; FIT coefficients 
  ;; for n=1 to n=4:
  ;;1=s, 2=p, 3=d, 4=f
  E=fltarr(8,4)
  E[*,0]= [-0.798, 1.14,  -1.24, -0.321,  0.173,  -0.0128,-0.0195, -0.0199  ]
  E[*,1]= [-2.78, -2.,    -1.09, -0.208, -0.127,  -0.136,  0.0375,  0.      ]
  E[*,2]= [-0.926,-0.344, -1.13, -0.164,  0.135,  -0.0675,-0.0372,  0.00486 ]
  E[*,3]= [-6.95, -3.49,  -0.617, 0.297,  0.174,   0.0395,-0.0364,  0.000221]
 
  if m le 2 then return,0.
  sigma=fltarr(4)
  for ii=0,3 do begin
     if m eq 3 then A=reform(D[*,ii])
     if m eq 4 then A=reform(E[*,ii])
     sigma[ii]=exp(A[0]/2. +total(A[1:4]*C[1:4]))*1.e-18
  endfor
  return, total(sigma)
end


@all_tables.pro
PRO neut_table,plot=plot,ps=ps
;; charge exchange n and m resolved charge exchange cross-sections
;; H(+)+H(n)-->H(m)+H(+)
;; The data for n=1 to n=3 is taken from ADAS. 
;; Higher states are derived by using the Principle of detailed balance
;; and by unbundeling the total cross sections from janev !!
  nmax=6        
  maxeb=800. 
  neb=1601


  ;; energy array
  ebarray=maxeb*dindgen(neb)/(neb-1.)
  deb=ebarray[2]-ebarray[1]
  print,'d-energy:', deb,'[kev]'
  ;; take cross sections from 'Mon. Not. R. Astron. Soc'
  ;; Excitation and charge transfer in H-H+ collisions at 5-80keV and
  ;; application to astrophysical shocks
  sigma_astro=fltarr(nmax,neb) 
  for ie=1L,neb-1 do begin
     for m=0,nmax-1 do begin
        sigma_astro[m,ie]=sigma1_astro(ebarray[ie],m+1)
     endfor
  endfor

  ;; cross sections from ADAS/JANEV
  mmax=nmax
  sigma_sum=fltarr(nmax,neb)   
  sigma=fltarr(mmax,nmax,neb)  
  ben=sigma_cx(12,1,1)
  print, ben
  for n=0,nmax-1 do begin
     sigma_sum[n,*]=sigma_cx(ebarray[*],n+1,-1)
     for m=0,nmax-1 do begin
        sigma[m,n,*]=sigma_cx(ebarray,n+1,m+1)
     endfor
  endfor



  ;; Plot the cross sections
  if keyword_set(plot) then begin
     loadct,39         ;;green=150;blue=50;yellow=200;red=254;black=0;white=255
     tvlct,0,150,0,150 ;; darker green!
     linthick=1.
     tthick=2.
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
     colorarr=[0,20,40,60,80,100,120,140,160,180,200,220,250]


     xran=[1.e-2,1.e3]
     yran=[1.e-20,1.e-12]
     xtit='Energy [keV/amu]'
     ytit='Cross-section [cm!u2!n]'
     for n=0L,nmax-1 do begin         ;; lower state
        file='PLOTS/neut_rate_n'+strtrim(string(n+1,f='(1i2)'),2)+'.eps'
        if keyword_set(ps) then device,filename=file
        ;if not(keyword_set(ps)) then window,n
        plot,[0.],/nodata,xran=xran,yran=yran,/ystyle,/ylog,/xlog $
             ,xtit=xtit,ytit=ytit,xthick=linthick,ythick=linthick,thick=linthick
        for i=0,12 do oplot,xran,[yran[0]*10.^i,yran[0]*10.^i] $
                            ,col=0,linesty=1,thick=1
        for i=0,5 do oplot,[xran[0]*10^i,xran[0]*10^i],yran $
                           ,col=0,linesty=1,thick=1
        xyouts,0.22,0.86,'Charge exchange from n='+string(n+1,f='(1i2)')+' to:',/norm
    
        for m=0,nmax-1 do begin      ;; upper state
         ;  if n eq 0 then oplot,ebarray[*],sigma_astro[m,*] $
         ;                       ,color=colorarr[m],thick=2
           oplot,ebarray[*],sigma[m,n,*],color=colorarr[m];,linesty=1
           print, m+1,n+1,max(sigma[m,n,*])
           if max(sigma[m,n,*]) gt 0. then begin
              
              xyouts,0.24,0.84-0.05*(m+1.) $
                     , 'm='+string(m+1,f='(1i2)') $
                     ,/norm,color=colorarr[m]
           endif
        endfor
        oplot,ebarray[*],sigma_sum[n,*],color=0,thick=tthick,linesty=1 
        dummy=''
        if not(keyword_set(ps)) then read,dummy
     endfor 
     if keyword_set(ps) then begin
        device, /close
        set_plot,'X'
        device,decomposed=0
     endif
     stop
  endif
  ;; write the cross sections into a file
  file='neuttable.bin'
  openw, lun, file, /get_lun
  writeu,lun, long(neb)
  writeu,lun, double(deb) 
  writeu,lun, long(nmax) 
  for ie=0L,neb-1 do begin 
     for n=0,nmax-1 do begin
        for m=0,nmax-1 do begin
              writeu, lun, double(sigma[m,n,ie])
        endfor
     endfor
  endfor
  close,lun
  free_lun, lun
  print, 'neutralization rate table written to:', file
end
