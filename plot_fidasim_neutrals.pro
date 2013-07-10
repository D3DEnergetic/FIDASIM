pro plot_fidasim_neutrals,ps=ps,only_halo=only_halo,only_beam=only_beam,loga=loga
  nlevs=15.
  linthick=1.
  set_plot,'X' & device, decomposed=0
  !P.charthick=1. & !P.charsize=1.5 & !P.thick=1.0
  if keyword_set(ps) then begin  
     set_plot, 'ps'

     device,  font_size=8, inches=0 , /encaps $
              ,xsize=7.9, ysize=6, /color, bits_per_pixel=8 
     device, /helvetica,font_index=3
     device, /symbol,font_index=4                                    
     device, /helvetica
     linthick=1.5
     leg_charsize=0.8
     !P.charsize=1. & !P.charthick=3. & !P.thick=1.5 & !P.font=0
  endif 
  loadct,39 ;; green=150;blue=50;yellow=200;red=254;black=0;white=255
  for i=0,64 do tvlct,255-i*3.3,255-i*3.3,255,i
  tvlct,0,0,0,0
  !P.background=255 & !P.color=0 &  !p.multi=0 


  path=dialog_pickfile(path='RESULTS/',/directory)

  runid=strsplit(path,'/',/extract,count=nid)
  print, runid
  runid=runid[nid-1]
  
  load_fidasim_results, fidasim,path
  shot=fidasim.inputs.shot


  xxc=fidasim.coords.xxc/100.
  yyc=fidasim.coords.yyc/100.
  zzc=fidasim.coords.zzc/100.
  xx=fidasim.coords.xx/100.
  yy=fidasim.coords.yy/100.
  zz=fidasim.coords.zz/100.
  dx=fidasim.coords.dx/100.
  dy=fidasim.coords.dy/100.
  dz=fidasim.coords.dz/100.
  nx=float(fidasim.coords.nx)
  ny=float(fidasim.coords.ny)
  nz=float(fidasim.coords.nz)
  xyzhead=fidasim.los.xyzhead/100.
  xyzlos=fidasim.los.xyzlos/100.
  rlos=sqrt(xyzlos[*,0]^2+xyzlos[*,1]^2)
  nchan=fidasim.los.nchan

  states=[0]
  halodens=total(fidasim.neutrals.halodens[*,*,*,*],4)
  fdens=total( fidasim.neutrals.fdens[*,*,*,*],4)
  hdens=total(fidasim.neutrals.hdens[*,*,*,*],4)
  tdens=total(fidasim.neutrals.tdens[*,*,*,*],4)
  nbidens= fdens+hdens+tdens
  if keyword_set(only_halo) then begin
     tot_dens=halodens
  endif else begin
     if keyword_set(only_beam) then begin
        tot_dens=nbidens
     endif else begin
        tot_dens= halodens+nbidens
     endelse
  endelse
  tot_dens_plasma=tot_dens
  ;; set density outside separatrix to 0.!
  for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
        for iz=0,nz-1 do begin
           if fidasim.plasma.rho[ix,iy,iz] gt 1. then  $
              tot_dens_plasma[ix,iy,iz]=0.
        endfor
     endfor
  endfor


  ;; PLOT TOP-DOWN VIEW
  file='PLOTS/top_down_beam_and_halo_neutrals.eps'
  if keyword_set(only_beam) then file='top_down_beam_neutrals.eps'
  if keyword_set(only_halo) then file='top_down_halo_neutrals.eps'
  if keyword_set(ps) then device,filename='PLOTS/'+file  $
  else window,1,xsize=700,ysize=600
  

  xran=[-2.5,2.5]
  yran=[-2.5,2.5]
  xticks=4
  yticks=4
  xtit='X [cm]'
  ytit='Y [cm]'
  ;; plot density profile
  plot,[0.], /nodata, xran=xran, yran=yran $
       ,xtit=xtit,ytit=ytit,/ystyle,/xstyle,xmar=xmar $
       ,ymar=ymar,/isotropic,xticks=xticks,yticks=yticks

  if keyword_set(transp) then xyouts,0.08,0.93,'TRANSP: Birth profile' $
                                     ,/norm,color=50
  pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
  
  tt='density of beam and halo neutrals [cm!u2!n]'
  if keyword_set(only_beam) then $
     tt='density of beam neutrals [cm!u2!n]'
  if keyword_set(only_halo) then $
     tt='density of halo neutrals [cm!u2!n]'
  xyouts,pos[2]+0.2,pos[3]+0.03,tt,/norm,alignm=1.
  ;; plot sectors and vessel
  i = 11
  x = [0.8,2.6]
  y = [-1.4,1.4]
  ri=1.1
  ra=2.2
  steps = 100
  phi = findgen(steps)/(steps-1)*2*!pi
  oplot,/polar,replicate(ri,steps),phi,thick=1
  oplot,/polar,replicate(ra,steps),phi,thick=1
  n_tor = 16
  density=reform(total(tot_dens,3))

     legpos=[.93,!y.window[0],.98,!y.window[1]]
     loadct,39
     for i=0,64 do tvlct,255-i*3.3,255-i*3.3,255,i
     tvlct,255,255,255,0

     if keyword_set(loga) then  density=alog10(density)>1.
     shade_surf,density,xxc,yyc,/noerase,  $
                ax=90,az=0,shades=bytscl(density,/NaN,top=253),$
                xtit=xtit,ytit=ytit,xran=xran,yran=yran,$ 
                zstyle=5,xmar=xmar,ymar=ymar,$
                /ystyle,/xstyle,pos=pos,xticks=xticks,yticks=yticks
     ;; plot colorbar
     mi = min(density,/NaN)
     ma = max(density,/NaN)
     xleg = [0,1]
     yleg = findgen(256)/255*(ma-mi)+mi
     dumarr = fltarr(2,256)
     dumarr(0,*) = yleg
     dumarr(1,*) = yleg
     legticks=5.
     ytickname = strarr(legticks)
     for i = 0,legticks-1 do begin
        yvalue = double(ma-mi)/(legticks-1.)*double(i)+mi

        yvalue=round(10.*yvalue/10.^fix(alog10(yvalue))) $
               *10.^fix(alog10(yvalue))/10.

        if keyword_set(loga) then yvalue=10.^yvalue
        ytickname(i) = string(yvalue,f='(1e7.1)')
     endfor
     shade_surf,dumarr,xleg,yleg,/noerase, $
                ax=90,az=0,shades=bytscl(dumarr,top=253,/NaN),$
                yran=[mi,ma],zstyle=5, /ystyle,$
                position=legpos
     tvlct,0,0,0,0
     plot,[0.],/nodata,xran=[0,1],yran=[mi,ma],/noerase  $
          ,pos=legpos, $
          xticks=1,yticks=legticks-1,charsize=leg_charsize,$
          xticklen=0.0001,yticklen=-0.01,$
          xtickname=(replicate(' ',2)),$
          ytickname=ytickname,xthick=linthick,ythick=linthick 
  




  plot,[0.], /nodata, xran=xran, yran=yran,/noerase $
       ,xtit=xtit,ytit=ytit,/ystyle,/xstyle,xmar=xmar $
       ,ymar=ymar,pos=pos,xticks=xticks,yticks=yticks $
       ,xthick=linthick,ythick=linthick
    tvlct,150,150,150,0
     for ii=0,nx-1 do begin
           oplot,[xx[ii],xx[ii]], [yy[0],yy[ny-1]+dy],thick=0.3*linthick,col=0
     endfor
     oplot,[xx[nx-1]+dx,xx[nx-1]+dx], [yy[0],yy[ny-1]+dy ],thick=0.3*linthick,col=0
     for ii=0,ny-1 do begin
           oplot,[xx[0],xx[nx-1]+dx], [yy[ii],yy[ii]],thick=0.3*linthick,col=0
     endfor
     oplot,[xx[0],xx[nx-1]+dx], [yy[ny-1]+dy,yy[ny-1]+dy],thick=0.3*linthick,col=0
     tvlct,0,0,0,0
 








  ;; POLOIDAL OVERVIEW
  file='PLOTS/neutrals_radial.eps'
  if keyword_set(only_beam) then file='radial_beam_neutrals.eps'
  if keyword_set(only_halo) then file='radial_halo_neutrals.eps'
  if keyword_set(ps) then device,filename='PLOTS/'+ file  $
  else window,3,xsize=700,ysize=600 
  ;; map density on a R-grid (rrc)
  rgrid=fltarr(nx,ny)
  for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
        rgrid[ix,iy]=sqrt(xxc[ix]^2+yyc[iy]^2)
     endfor
  endfor
  nr=nx
  rrc=dindgen(nr)/(nr-1.)*(max(rgrid)-min(rgrid))+min(rgrid)
  dr=rrc[1]-rrc[0]
  density=replicate(0.,nr,nz)
  density_plasma=replicate(0.,nr,nz)
  for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
        dummy=min(abs(rrc-rgrid[ix,iy]),ir)
        for iz=0,nz-1 do begin
           density[ir,iz]=density[ir,iz]+tot_dens[ix,iy,iz]
           density_plasma[ir,iz]=density_plasma[ir,iz]+tot_dens_plasma[ix,iy,iz]
        endfor
     endfor
  endfor
   
  ztit='Z [m]'
  rtit='R [m]'
  xmar=[8.0,9.]
  rran=[.8,2.5]
  zran=[-1.,1.]

  title='XY-integrated density [1/cm^2]' 
  tvlct,0,0,0,0
  plot,[0.], /nodata, xran=rran, yran=zran $
       ,xtit=rtit, ytit=ztit,/ystyle,/xstyle,/isotropic $
       ,xticks=xticks,yticks=yticks,xmar=xmar
  index=where(rrc le rran[1],nind)
  if nind gt 0 then begin
     density=density[index,*]
     rrc=rrc[index]
     nr=nind
  endif
  help, density
  pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
  xyouts,pos[2]+0.2,0.94,'density [cm!u2!n]',/norm,alignm=1.

  legpos=[pos[2]+0.15,!y.window[0],pos[2]+0.2,!y.window[1]]
  loadct,39
  for i=0,64 do tvlct,255-i*3.3,255-i*3.3,255,i
  tvlct,255,255,255,0
  if keyword_set(loga) then  density=alog10(density)>1.
  
  shade_surf,density,rrc,zzc,/noerase $
             ,ax=90,az=0,shades=bytscl(density,/NaN,top=253) $
             ,xtit=rtit,ytit=ztit,xran=rran,yran=zran $ 
             ,zstyle=5,/ystyle,/xstyle,pos=pos $
             ,xticks=xticks,yticks=yticks
  mi = min(density,/NaN)
  ma = max(density,/NaN)
  xleg = [0,1]
  yleg = findgen(256)/255*(ma-mi)+mi
  dumarr = fltarr(2,256)
  dumarr(0,*) = yleg
  dumarr(1,*) = yleg
  legticks=5.
  ytickname = strarr(legticks)
  for i = 0,legticks-1 do begin
     yvalue = double(ma-mi)/(legticks-1.)*double(i)+mi
     yvalue=round(10.*yvalue/10.^fix(alog10(yvalue))) $
            *10.^fix(alog10(yvalue))/10.
     if keyword_set(loga) then yvalue=10.^yvalue
     ytickname(i) = string(yvalue,f='(1e7.1)')
  endfor
  shade_surf,dumarr,xleg,yleg,/noerase, $
             ax=90,az=0,shades=bytscl(dumarr,top=253,/NaN),$
             yran=[mi,ma],zstyle=5, /ystyle,$
             position=legpos,charsize=leg_charsize
  tvlct,0,0,0,0
  plot,[0.],/nodata,xran=[0,1],yran=[mi,ma],/noerase  $
       ,pos=legpos, $
       xticks=1,yticks=legticks-1,charsize=leg_charsize, $
       xticklen=0.0001,yticklen=-0.01,$
       xtickname=(replicate(' ',2)),$
       ytickname=ytickname,xthick=linthick,ythick=linthick
  
  
  plot,[0.], /nodata, xran=rran, yran=zran,/noerase $
       ,xtit=rtit, ytit=ztit,/ystyle,/xstyle,pos=pos $
       ,xticks=xticks,yticks=yticks,xthick=linthick,ythick=linthick

  
  
  tvlct,155,155,155,0
  rr=rrc-0.5*dr
  for ii=0,nr-1 do begin
     oplot,[rr[ii],rr[ii]], [zz[0],zz[nz-1]+dz],thick=0.5*linthick,col=0
  endfor
  oplot,[rr[nr-1]+dr,rr[nr-1]+dr], [zz[0],zz[nz-1]+dz],thick=0.5*linthick,col=0
  for ii=0,nz-1 do begin
     oplot,[rr[0],rr[nr-1]+dr], [zz[ii],zz[ii]],thick=0.5*linthick,col=0
  endfor
  oplot,[rr[0],rr[nr-1]+dr], [zz[nz-1]+dz,zz[nz-1]+dz],thick=0.5*linthick,col=0
  tvlct,0,0,0,0
  
  
  
  if keyword_set(ps) then device, /close
  stop
end                             ;of programm

