pro rotate_uvw,uvw,Arot,Brot,Crot,updown,xyz
  A2rot=reform(Arot[0,*,*])
  B2rot=reform(Brot[0,*,*])
  C2rot=reform(Crot[0,*,*])
  if updown lt 0. then  qrz=MATRIX_MULTIPLY(A2rot,uvw) 
  if updown gt 0. then  qrz=MATRIX_MULTIPLY(B2rot,uvw)
  xyz=MATRIX_MULTIPLY(C2rot,qrz)
end


pro plot_dalpha_inputs,inputs,nbgeom,coords,plasma,detector,dens=dens
;  restore,'plot_dalpha_inputs.idl'
  simulate_nbi=0
  plot_grid=1
 ;;----------------------------------
 ;;;plot profiles
 ;;----------------------------------
  tvlct,0,150,0,150
  if inputs.ps eq 1 then begin
     !P.multi=[0,0,0]
     set_plot,'ps'
     device,  font_size=8, inches=0 , /encaps $
                ,xsize=7.9, ysize=6, /color, bits_per_pixel=8 
     device, /helvetica,font_index=3
     device, /symbol,font_index=4                                    
     device, /helvetica
     linthick=1.5
     !P.charsize=1. & !P.charthick=3. & !P.thick=linthick & !P.font=0
     file='PLOTS/input_profiles_'+string(inputs.shot,f='(i5)')+'.eps'
     device,filename=file;,xsize=15.8,ysize=12
  endif else begin
     set_plot,'X'
     device,decomposed=0
     window,1,xsize=790,ysize=600
     linthick=1
  endelse
  xmar=[4.5,.5]
  ymar=[1.5,.5]
  xticks=5
  !P.background=255 & !P.color=0 & !P.multi=[0,2,2]
 

  index=long(indgen(400)/399.*coords.ng)
  rho=plasma.rho_grid[index]
  sind=sort(rho)
  rho=rho[sind]
  index=index[sind]
  ;;PLOT TE and TI
  xtit='rho poloidal'
  xran=[0,1.0]
  plot,rho,plasma.te[index],xran=xran,/xsty,yran=[0.,6],/ysty $
       ,ytit='temperature [keV]',xmargin=xmar,ymargin=ymar,xthick=linthick $
       ,ythick=linthick,xticks=xticks,xtickname=(replicate(' ',xticks+1)),yticks=3
  oplot,rho,plasma.ti[index],color=254
  xyouts,0.25,0.85,'T!de!n,',/norm, color=1
  xyouts,0.3,0.85,'T!di!n',/norm, color=254
  xyouts,0.745,0.93,string(inputs.shot,'(i5)')+'@' $
         +string(inputs.time,'(1f5.3)')+'s',/norm
  ;; PLOT the DENSITIES
  dene=plasma.dene[index]*1.d6/1.e19
  plot,rho,dene,yran=[0.,8],xran=xran,/xsty,/ysty $
       ,xthick=linthick,ythick=linthick $
       ,ytit='density [10!u19!nm!u-!u3!n]' $
       ,xmargin=xmar,ymargin=ymar,xticks=xticks $
       ,xtickname=(replicate(' ',xticks+1)),yticks=4
  oplot,rho,plasma.deni[index]*10.d6/1.e19,color=50
  oplot,rho,plasma.denp[index]*1.d6/1.e19,color=200


  index2=where(coords.r_grid gt 165. and coords.zc gt -10. and coords.zc lt 10.)
  oplot,plasma.rho_grid[index2],plasma.denf[index2]*10.d6/1.e19,color=254,psym=1,symsize=0.5
  xyouts,0.65,0.85 ,'n!de!n',/norm
  xyouts,0.7,0.85,'n!dD!n',/norm,color=200
  xyouts,0.75,0.85 ,'n!dimp!n*10',/norm,color=50
  xyouts,0.85,0.85,'n!dfi!n*10',/norm,color=254 
 ; xyouts,0.57,0.95,'x 10!u19!n',/norm,color=0,chars=.9
  ;; PLOT ZEFF
  ymar=[3.,-1.]
  plot,[0.],/nodata,yran=[1.,2.],/ysty,xran=xran,/xsty,xthick=linthick $
       ,ythick=linthick,yticks=2 $
       ,xtit=xtit,xmargin=xmar,ymargin=ymar,xticks=xticks
  oplot,rho,plasma.zeff[index]
  xyouts,0.13,0.35,'effective charge Z!deff!n',/norm,color=0
  ;; PLOT VTOR
  plot,[0.],/nodata,yran=[0.d0,250],/ysty,xran=xran,/xsty,xthick=linthick $
       ,ythick=linthick,yticks=2 $
       ,xtit=xtit,xmargin=xmar,ymargin=ymar,xticks=xticks
  oplot,rho,plasma.vtor[index]/1.e5 
  xyouts,0.61,0.35,'toroidal rotation [km/s]',/norm,color=0
  ymar=[3.0,2.]

  ;; ----------------------------------
  ;; plot geometry overview
  ;; ----------------------------------
  ;; calculated parameters
  xx=coords.xx/100. & yy=coords.yy/100. & zz=coords.zz/100.
  dx=coords.dx/100. & dy=coords.dy/100. & dz=coords.dz/100.
  nx=coords.nx & ny=coords.ny & nz=coords.nz
  coords.xc=coords.xc/100.
  coords.yc=coords.yc/100.
  coords.zc=coords.zc/100.
  coords.x=coords.x/100.
  coords.y=coords.y/100.
  coords.z=coords.z/100.
  coords.xxc=coords.xxc/100.
  coords.yyc=coords.yyc/100.
  coords.zzc=coords.zzc/100.
  wz=where(abs(coords.z) eq min(abs(zz)) and plasma.rho_grid ge 1)
  wy=where(abs(coords.y) eq min(abs(yy)) and plasma.rho_grid ge 1.)

  ;; PLOT SETTINGS
  if inputs.rotate eq 0 then begin
     if inputs.isource le 3 then begin
        xran=[-0.,2.5]
        yran=[-3.,.5]
        if inputs.npa eq 1 then begin 
           xran=[0,4.]
           yran=[-2.5,2.5]
        endif
        if simulate_nbi eq 1 then xran=[-3.,7.]
        if simulate_nbi eq 1 then yran=[-7.,3.]
     endif else begin
        xran=[-2.5,0]
        yran=[-.5,2.5]
     endelse 
  endif else begin
     xran=[-1,2.5]
     yran=[-1.,2.5]
  endelse
  if inputs.isource gt 3 then begin
  endif

  ;; TOP DOWN OVERVIEW
  file='PLOTS/input_top_down_view.eps'
  if inputs.ps then device,filename=file,xsize=7.9, ysize=6
  !p.multi=[0,0,0]
  if not(inputs.ps) then window,2,xsize=500,ysize=500
  plot,[0.], /nodata, xran=xran, yran=yran,/isotropic $
       ,xtit='X [m]',ytit='Y [m]',/ystyle,/xstyle $
       ,xthick=linthick,ythick=linthick,thick=linthick $
       ,xmar=xmar, ymar=ymar
 ;; Fast-ion distribution
  if total(plasma.denf)gt 0. then begin
     pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
     for i=0,20 do tvlct,255-i*10,255-i*10,255-i*10,i
     
     f3d=fltarr(coords.nx,coords.ny,coords.nz)*0.d0
     for i=0L,coords.nx-1 do begin
        for j=0,coords.ny-1 do begin
           for k=0L,coords.nz-1 do begin
              jj=i+coords.nx*j+coords.nx*coords.ny*k
              f3d[i,j,k]=plasma.denf[jj]
           endfor
        endfor
     endfor 
     index=where(coords.zzc gt -0.2 and coords.zzc lt .2)
     for i=0,60 do tvlct,255-i*4.,255-i*4.,255,i
     tvlct,0,0,0,254
     ben=reform(total(f3d[*,*,index],3))
     shade_surf,ben,coords.xxc,coords.yyc,/noerase, $ 
                ax=90,az=0,shades=bytscl(ben,top=253,/NaN),$
                yran=yran,xran=xran, $ 
                zstyle=5,/ystyle,/xstyle, $
                position=pos
     tvlct,0,0,0,0
     plot,[0.], /nodata, xran=xran, yran=yran $
       ,xtit='X [m]',ytit='Y [m]',/ystyle,/xstyle $
       ,xthick=linthick,ythick=linthick,thick=linthick $
       ,xmar=xmar, ymar=ymar,/noerase,pos=pos
  endif
;; Plot grid!
  ;loadct,39
  tvlct,150,150,150,0
  if simulate_nbi eq 1 then begin
     for ii=0,nx-1 do begin
        if ii mod 5 eq 0 then begin
           oplot,[xx[ii],xx[ii]], [yy[0],yy[ny-1]+dy],thick=0.5*linthick
        endif
     endfor
     oplot,[xx[nx-1]+dx,xx[nx-1]+dx], [yy[0],yy[ny-1]+dy ],thick=0.5*linthick
     for ii=0,ny-1 do begin
        if ii mod 5 eq 0 then begin
           oplot,[xx[0],xx[nx-1]+dx], [yy[ii],yy[ii]],thick=0.5*linthick
        endif
     endfor
     oplot,[xx[0],xx[nx-1]+dx], [yy[ny-1]+dy,yy[ny-1]+dy],thick=0.5*linthick
  endif else begin
     if plot_grid eq 1 then begin
     for ii=0,nx-1 do oplot,[xx[ii],xx[ii]], [yy[0],yy[ny-1]+dy],color=0,thick=0.5*linthick
     oplot,[xx[nx-1]+dx,xx[nx-1]+dx], [yy[0],yy[ny-1]+dy ],color=0,thick=0.5*linthick
     for ii=0,ny-1 do oplot,[xx[0],xx[nx-1]+dx], [yy[ii],yy[ii]],color=0,thick=0.5*linthick
     oplot,[xx[0],xx[nx-1]+dx], [yy[ny-1]+dy,yy[ny-1]+dy],color=0,thick=0.5*linthick
     endif
   endelse
     tvlct,0,0,0,0
  i = 11
  x = [0.8,2.6]
  y = [-1.4,1.4]
  ri=1.1
  ra=2.2
  steps = 100
  phi = findgen(steps)/(steps-1)*2*!pi
  oplot,/polar,replicate(ri,steps),phi, thick=1.2
  oplot,/polar,replicate(ra,steps),phi, thick=1.2


     ;; calucalte NBI ray
     isource=coords.isource
     xyz_start_src=nbgeom.xyz_src[isource,*]/100.
     du=-sqrt(xyz_start_src[0]^2+xyz_start_src[1]^2+xyz_start_src[2]^2)
     uvw_ray=[du,0.d0,0.d0]
     updown=-1
     rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*],updown, xyz_ray
     xyz_end_src=nbgeom.xyz_src[isource,*]/100.  + xyz_ray
     oplot,[xyz_start_src[0],xyz_end_src[0]] $
        ,[xyz_start_src[1],xyz_end_src[1]],thic=4,color=0
     



     

     for chan=0, detector.det.nchan-1 do begin
        if detector.det.nchan lt 100 or (chan mod 10) eq 0 then begin
           if inputs.calc_wght[0] eq 1 and inputs.ichan_wght gt 0 then begin
              if chan ne inputs.ichan_wght-1 then continue
           endif

           losv=[detector.det.xhead[chan]-detector.det.xlos[chan], $
                    detector.det.yhead[chan]-detector.det.ylos[chan], $
                    detector.det.zhead[chan]-detector.det.zlos[chan]]/100.
           
           an=acos((xyz_ray[0]*losv[0]+xyz_ray[1]*losv[1]+xyz_ray[2]*losv[2]) $
                 /sqrt(xyz_ray[0]^2 + xyz_ray[1]^2 + xyz_ray[2]^2) $
                 /sqrt(losv[0]^2 + losv[1]^2 + losv[2]^2))
           print, 'angle between LOS and beam: ' , an*!radeg
           xhead=double(detector.det.xhead[chan])/100.d0
           yhead=double(detector.det.yhead[chan])/100.d0
           zhead=double(detector.det.zhead[chan])/100.d0
           xlos=double(detector.det.xlos[chan])/100.d0
           ylos=double(detector.det.ylos[chan])/100.d0
           if (yhead-ylos) lt 0. then dyy=7.d0 else dyy=-7.d0
           dxx=dyy * double(xlos-xhead)/double(ylos-yhead)
           if zhead lt 0.5 then col=150 else col=50
           oplot,[xhead,xhead+dxx],[yhead,yhead+dyy],color=col,thick=1.
        endif
     endfor  



  ;; POLOIDAL OVERVIEW
  xran=[.8,2.5]
  if simulate_nbi eq 1 then xran=[0.8,9.5]
  zran=[-1.5,1.5]
  file='PLOTS/input_radial_view.eps'
  if inputs.ps then device,filename=file
  !p.multi=[0,0,0]
  if not(inputs.ps) then window,3,xsize=900,ysize=500
  plot,[0.], /nodata, xran=xran, yran=zran,/isotropic $
       ,xtit='R [m]', ytit='Z [m]',/ystyle,/xstyle $
       ,xthick=linthick,ythick=linthick,thick=linthick $
       ,xmar=xmar, ymar=ymar


  if total(plasma.denf)gt 0. then begin
     pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
     nr=coords.nx
     rrc=(dindgen(nr)/(nr-1.)*(max(coords.r_grid)-min(coords.r_grid)) $
          +min(coords.r_grid))/100.
     ben=fltarr(nr,coords.nz)*0.d0
     for i=0L,coords.nx-1 do begin
        for j=0,coords.ny-1 do begin
           for k=0,coords.nz-1 do begin
              dummy=min(abs(sqrt(coords.xxc[i]^2+coords.yyc[j]^2)-rrc),index)
              ben[index,k]=f3d[i,j,k]
           endfor
        endfor 
     endfor
     tvlct,254,254,254,0
     shade_surf,ben,rrc,coords.zzc,/noerase, $ 
                ax=90,az=0,shades=bytscl(ben,top=253,/NaN),$
                yran=zran,xran=xran, $ 
                zstyle=5,/ystyle,/xstyle, $
                position=pos
     tvlct,0,0,0,0
     plot,[0.], /nodata, xran=xran, yran=zran $
          ,xtit='R [m]', ytit='Z [m]',/ystyle,/xstyle $
          ,xthick=linthick,ythick=linthick,thick=linthick $
          ,xmar=xmar, ymar=ymar,/noerase,pos=pos
  endif


 
     tvlct,150,150,150,0
     rrmin=min(coords.r_grid)/100.
     rrmax=max(coords.r_grid)/100.
     nr=nx
     rrc=dindgen(nr)/(nr-1.)*(rrmax-rrmin)+rrmin
     dr=rrc[1]-rrc[0]
     rr=rrc-0.5*dr

     rr=rrc-0.5*dr
     if simulate_nbi eq 1 then begin
        for ii=0,nr-1 do begin
           if ii mod 5 eq 0 then begin
              oplot,[rr[ii],rr[ii]], [zz[0],zz[nz-1]+dz],thick=0.5*linthick
           endif
        endfor
        oplot,[rr[nr-1]+dr,rr[nr-1]+dr], [zz[0],zz[nz-1]+dz],thick=0.5*linthick
        for ii=0,nz-1 do begin
           if ii mod 4 eq 0 then begin
              oplot,[rr[0],rr[nr-1]+dr], [zz[ii],zz[ii]],thick=0.5*linthick
           endif
        endfor
        oplot,[rr[0],rr[nr-1]+dr], [zz[nz-1]+dz,zz[nz-1]+dz],thick=0.5*linthick
     endif else begin
        if plot_grid eq 1 then begin
        for ii=0,nr-1 do oplot,[rr[ii],rr[ii]], [zz[0],zz[nz-1]+dz],thick=0.5*linthick
        oplot,[rr[nr-1]+dr,rr[nr-1]+dr], [zz[0],zz[nz-1]+dz],thick=0.5*linthick
        for ii=0,nz-1 do oplot,[rr[0],rr[nr-1]+dr], [zz[ii],zz[ii]],thick=0.5*linthick
        oplot,[rr[0],rr[nr-1]+dr], [zz[nz-1]+dz,zz[nz-1]+dz],thick=0.5*linthick
        endif
     endelse
     tvlct,0,0,0,0  


  ;; plot NBI
  ;if not keyword_set(dens) then begin ;; no passive fida simulation
 
  isource=coords.isource
  xyz_start_src=nbgeom.xyz_src[isource,*]/100.
     nii=400
     for ii=0,nii-1 do begin
        uvw_ray=[-1.d0,0.d0,0.d0]*(5.+float(ii)/float(nii-1)*4.)
        rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                   ,nbgeom.Brot[isource,*,*] $
                   ,nbgeom.Crot[isource,*,*] $
                   ,-1,nbi_ray 
        xyz_end_src=nbgeom.xyz_src[isource,*]/100.  + nbi_ray
        r_end=sqrt(xyz_end_src[0]^2+xyz_end_src[1]^2)
        z_end=xyz_end_src[2]
        if ii ge 1 then oplot,[r_endb1,r_end] , [z_endb1,z_end],thic=4,color=0
        r_endb1=r_end
        z_endb1=z_end
        rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                   ,nbgeom.Brot[isource,*,*] $
                   ,nbgeom.Crot[isource,*,*] $
                   ,1,nbi_ray
        xyz_end_src=nbgeom.xyz_src[isource,*]/100.  + nbi_ray
        r_end=sqrt(xyz_end_src[0]^2+xyz_end_src[1]^2)
        z_end=xyz_end_src[2]
        if ii ge 1 then oplot,[r_endb2,r_end] , [z_endb2,z_end],thic=4,color=0
        r_endb2=r_end
        z_endb2=z_end
     endfor

     for chan=0, detector.det.nchan-1 do begin
        if detector.det.nchan lt 100 or (chan mod 10) eq 0 then begin
           if inputs.calc_wght[0] eq 1 and inputs.ichan_wght gt 0 then begin
              if chan ne inputs.ichan_wght-1 then continue
           endif
           xhead=detector.det.xhead[chan]/100.
           yhead=detector.det.yhead[chan]/100.
           zhead=detector.det.zhead[chan]/100.
           rhead=sqrt(xhead^2+yhead^2)
           xlos=detector.det.xlos[chan]/100.
           ylos=detector.det.ylos[chan]/100.
           zlos=detector.det.zlos[chan]/100.
           dx=xhead - xlos
           dy=yhead - ylos
           dz=zhead - zlos
           rb=rhead
           zb=zhead
           for ii=1.,10 do begin
              if zhead gt 1. then begin
                 dii=-1*ii/5.
                 x=xhead+dx/dz*dii
                 y=yhead+dy/dz*dii
                 z=zhead+dz/dz*dii
                 r=sqrt(x^2+y^2)
                 oplot,[rb,r],[zb,z],color=ii*20,thick=1.
                 rb=r
                 zb=z
              endif else begin
                 dii=1*ii/5.
                 if (yhead-ylos) gt 0. then dii=-dii
                 x=xhead+dx/dy*dii
                 y=yhead+dy/dy*dii
                 z=zhead+dz/dy*dii
                 r=sqrt(x^2+y^2)
                 oplot,[rb,r],[zb,z],color=ii*20,thick=1.
                 rb=r
                 zb=z
              endelse
           endfor
        endif
     endfor
 
  if inputs.ps then begin
     file='PLOTS/inputs.eps'
     device,filename=file
  endif


end
