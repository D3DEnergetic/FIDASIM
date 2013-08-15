;;--------------------------------------------------
;;--------------- PLOT DALPHA INPUTS ---------------
;;--------------------------------------------------
;; ROUTINE of FIDASIM to plot the input kinetic profiles and the
;; applied geometry (simulation grid and LOS) before starting the
;; FORTRAN part.
;; written by Benedikt Geiger 2013
;; additional subroutines used: 
pro plot_dalpha_inputs,inputs,nbgeom,coords,plasma,detector,dens=dens
  plot_grid=1
  ab=1
  color_gradient=1
  tvlct,0,150,0,150
  if inputs.ps eq 1 then begin
     set_plot,'ps'
     device,font_size=8, inches=0,/encaps $
           ,xsize=8.5,ysize=5,/color,bits_per_pixel=8,/helvetica
     !P.charsize=1. & !P.charthick=1. & !P.thick=1.746 & !P.font=0
     file='PLOTS/fidasim_input_profiles_'+string(inputs.shot,f='(i5)')+'.eps'
     device,filename=inputs.root_dir+file
  endif else begin
     set_plot,'X'
     device,decomposed=0
     !P.charsize=2.5 & !P.charthick=1. & !P.thick=1. & !P.font=-1 
     window,2,xsize=850,ysize=500
  endelse
  chars=0.9*!p.charsize
  los_thick=2.
  !P.background=255 & !P.color=0 & !P.multi=[0,2,2]







  ;;----------------------------------
  ;;--------- plot profiles ----------
  ;;----------------------------------
  xtit='rho'
  xran=[0,1.1]
  ;; define grid cells whose profiles are plotted
  index=long(indgen(400)/399.*coords.ng)
  rho=plasma.rho_grid[index]
  sind=sort(rho)
  rho=rho[sind]
  index=index[sind]
  ;; define plot position:
  dxpos=0.4 ;; width of plot
  dx2r=0.02  ;; right margin 
  dypos=0.41 ;; height of plot
  dyu=0.02  ;; upper margin
  dx1r=(1.-2*dxpos-dx2r)/2.
  pos1=[1.-2*dxpos-dx1r-dx2r,1.-dypos-dyu    ,1.-dxpos-dx1r-dx2r,1.-dyu    ]
  pos2=[1.-dxpos-dx2r       ,pos1[1]         ,1.-dx2r       ,pos1[3]       ]
  pos3=[pos1[0]             ,1.-2*dypos-2*dyu,pos1[2]       ,1.-dypos-2*dyu]
  pos4=[pos2[0]             ,pos3[1]         ,pos2[2]       ,pos3[3]       ]
  ;; -------------------------------
  ;; ----- PLOT TEMPERATURE --------
  ;; -------------------------------
  plot,rho,plasma.te[index],pos=pos1 $
       ,xran=xran,/xsty $
       ,/xminor,/yminor $
       ,yran=[0,ceil(max(plasma.te))],/ysty,yticks=ceil(max(plasma.te)) $
       ,xtickname=(replicate(' ',20))
  ytit='keV'
  xyouts,0.04,(pos1[1]+pos1[3])/2.,ytit,/norm,orientation=90,align=0.5
  oplot,rho,plasma.ti[index],color=254
  xyouts,pos1[0]+0.20,pos1[3]-0.10,'T!de!n',/norm
  xyouts,pos1[0]+0.25,pos1[3]-0.10,'T!di!n',/norm,color=254
  ;; --------------------------
  ;; ----- PLOT DENSITY -------
  ;; --------------------------
  dene=plasma.dene[index]*1.d6/1.e19
  plot,rho,dene,pos=pos2,/noerase $
       ,xran=xran,/xsty $
       ,/xminor,/yminor $
       ,yran=[0.,ceil(max(dene))],/ysty,yticks=ceil(max(dene)) $
       ,xtickname=(replicate(' ',20))
  ytit='10!u19!nm!u-!u3!n' 
  xyouts,pos1[2]+0.04,(pos2[1]+pos2[3])/2.,ytit,/norm,orientation=90,align=0.5
  oplot,rho,plasma.deni[index]*1.d6/1.e19,color=50
  oplot,rho,plasma.denp[index]*1.d6/1.e19,color=200
  xyouts,pos2[0]+0.15,pos2[3]-0.1,'n!de!n',/norm,chars=chars
  xyouts,pos2[0]+0.20,pos2[3]-0.1,'n!dD!n',/norm,color=200,chars=chars
  xyouts,pos2[0]+0.25,pos2[3]-0.1,'n!dB!n',/norm,color=50 ,chars=chars
  if total(plasma.denf) gt 0 then begin
     index2=where(coords.rrc_grid gt 165.      $
                  and coords.zzc_grid gt -10.  $
                  and coords.zzc_grid lt  10.)
     oplot,plasma.rho_grid[index2],plasma.denf[index2]*1.d6/1.e19 $
           ,color=254,psym=1,symsize=0.5
     xyouts,pos2[0]+0.30,pos2[3]-0.1,'n!dfi!n',/norm,color=254,chars=chars
  endif

  ;; -----------------------
  ;; --- PLOT ZEFF ---------
  ;; -----------------------
  plot,rho,plasma.zeff[index],/noerase,pos=pos3 $
       ,xran=xran,/xsty $
       ,/xminor,/yminor $
       ,yran=[1.,ceil(max(plasma.zeff[index]))],/ysty,yticks=2
  ytit='Z!deff!n'
  xyouts,0.03,(pos3[1]+pos3[3])/2.,ytit,/norm,orientation=90,align=0.5
  xyouts,(pos3[0]+pos3[2])/2.,0.015,xtit,/norm,align=0.5 
  xyouts,0.5*(pos3[0]+pos3[2]),pos3[3]-0.10,'effective charge' $
         ,/norm,align=0.5,chars=chars
  ;; -----------------------
  ;; --- PLOT VTOR ---------
  ;; -----------------------
  plot,rho,plasma.vtor[index]/1.e5,/noerase,pos=pos4 $
       ,xran=xran,/xsty $
       ,/xminor,/yminor $
       ,yticks=2
  ytit='km/s'
  xyouts,pos1[2]+0.03,(pos4[1]+pos4[3])/2.,ytit,/norm,orientation=90,align=0.5
  xyouts,(pos4[0]+pos4[2])/2.,0.015,xtit,/norm,align=0.5 
  xyouts,0.5*(pos4[0]+pos4[2]),pos4[3]-0.1,'toroidal rotation' $
         ,/norm,align=0.5,chars=chars
  ;; ----------- print shot-number and time -------
  xyouts,pos2[2]-0.01,pos2[3]-0.05,string(inputs.shot,'(i5)')+'@' $
         +string(inputs.time,'(1f5.3)')+'s',/norm,/align,chars=chars






  ;; ---------------------------------------
  ;; ------------ PLOT GEOMETRY ------------
  ;; ---------------------------------------
  xx=coords.xx/100. & yy=coords.yy/100. & zz=coords.zz/100.
  dx=coords.dx/100. & dy=coords.dy/100. & dz=coords.dz/100.
  nx=coords.nx & ny=coords.ny & nz=coords.nz
  xxc=coords.xxc/100.
  yyc=coords.yyc/100.
  zzc=coords.zzc/100.
  ;; range settings
  if inputs.rotate eq 0 then begin
     if inputs.isource le 3 then begin
        xran=[-0.5,3.]
        yran=[-3.,1.]
        if inputs.npa eq 1 then begin 
           xran=[0,4.]
           yran=[-2.5,2.5]
        endif
     endif else begin
        xran=[-2.5,0]
        yran=[-.5,2.5]
     endelse 
  endif else begin
     xran=[-1,2.5]
     yran=[-1.,2.5]
  endelse
  rran=[.8,2.5]
  zran=[-1.2,1.2]
  ;; dummy plot to define plot positions
  if inputs.ps then device,filename=inputs.root_dir+'idl.eps',xsiz=8.5,ysiz=5. $
  else window,3,xsize=850,ysize=500
  !p.multi=[0,2,0]
  ymar=[3.0,1.5]
  xmar=[5.4,0]
  plot,[0.],/nodata,xran=xran,yran=yran,/isotropic $
      ,/ystyle,/xstyle $
      ,xmar=xmar,ymar=ymar
  dxpos=-0.015
  pos1=[!x.window[0]+dxpos,!y.window[0],!x.window[1]+dxpos,!y.window[1]]
  xmar=[5.9,2.1]
  plot,[0.], /nodata, xran=rran, yran=zran,/isotropic $
       ,/ystyle,/xstyle  $
       ,xmar=xmar,ymar=ymar
  pos2=[!x.window[0]+dxpos,!y.window[0],!x.window[1]+dxpos,!y.window[1]]
  file='PLOTS/fidasim_geom_'+string(inputs.shot,f='(i5)')+'.eps'
  if inputs.ps then device,filename=inputs.root_dir+file,xsize=8.5, ysize=5.


  ;; --------------------------------------
  ;; ----------- TOP DOWN VIEW ------------
  ;; --------------------------------------
  if total(plasma.denf)eq 0. then begin
     plot,[0.],/nodata,xran=xran, yran=yran $
             ,/ystyle,/xstyle $
             ,pos=pos1,yticks=4,xticks=3
  endif else begin
     for i=1,254 do tvlct,256-i,256-i,256-i,i
  ;; PLOT Fast-ion density
     index=where(zzc gt 0. and zzc lt .2,nind)
     zplot=reform(total(plasma.denf[*,*,index],3)/nind)*1.d6
     nlevels=8
     ma = max(plasma.denf,/NaN)*1.d6
     levelvals=(dindgen(nlevels)+0.0 )/nlevels*ma
     c_colors=dindgen(nlevels)*253./(nlevels-1.)+1
     contour,zplot,xxc,yyc, xran=xran, yran=yran $
             ,/ystyle,/xstyle $
             ,c_colors=c_colors        $
             ,levels=levelvals,/fill,pos=pos1,yticks=4,xticks=3
     loadct,39
  endelse
  xyouts,pos1[0]-0.053,(pos1[1]+pos1[3])/2.,'Y [m]' $
         ,/norm,orientation=90,align=0.5
  xyouts,(pos1[0]+pos1[2])/2.,0.02,'X [m]',/norm,align=0.5
  if plot_grid then oplot_fidasim_grid,coords,/tor,thick=0.5*!p.thick
 


  if total(plasma.denf)ne 0 then xyouts,pos1[2]-.007,pos1[3]-.05,'z=0.1m'+'!n' $
                                        ,/norm,chars=chars,align=1

  ;; plot magnetic field
  if 0 then begin
     counter=0L
     dummy=min(abs(zzc),zc)
     ind=where(plasma.rho_grid[*,*,zc] le 1.,nind)
     for i=0,nx-1 do begin
        for j=0,ny-1 do begin
           if plasma.rho_grid[i,j,zc] le 1. then begin
              if counter gt 20 then begin ;; do not oplot b for all cells
                 oplot,[xxc[i],xxc[i]+plasma.bx[i,j,zc]*.05]  $
                       ,[yyc[j],yyc[j]+plasma.by[i,j,zc]*.05] $
                       ,thick=2
                 oplot,[xxc[i]],[yyc[j]],psym=1,color=254
                 counter=0L
              endif
              counter++
           endif
        endfor
     endfor
  endif
  ;; PLOT NBI ray
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
        ,[xyz_start_src[1],xyz_end_src[1]],thic=4,color=200
  if isource lt 4 and inputs.rotate eq 0 then begin
     xyouts,1.9,-2.8,'NBI Q'+string(isource+1,f='(i1)'),color=200,chars=chars
  endif
 
  ;; PLOT LOS
  ;; CALCULATE THE LENGHT (dl*nstep) OF THE LOS TO BE DISPLAYED
  dl=10.   ;; [cm]
  ;; determine the index of the toroidal LOS
  index=where(detector.det.zhead lt 100. and detector.det.zhead gt -100.,ntor)
  if ntor gt 0 then begin
     dx=mean(detector.det.xhead[index]) - mean(detector.det.xlos[index])
     dy=mean(detector.det.yhead[index]) - mean(detector.det.ylos[index])
     dz=mean(detector.det.zhead[index]) - mean(detector.det.zlos[index])
     nstep_tor=sqrt(dx^2+dy^2+dz^2)/dl
     nstep2_tor=float(round(nstep_tor*1.2))
  endif
  ;; determine the index of the poloidal LOS (viewing from top to bottom)
  index=where(abs(detector.det.zhead) gt 100.,ntor)
  if ntor gt 0 then begin
     dx=mean(detector.det.xhead[index]) - mean(detector.det.xlos[index])
     dy=mean(detector.det.yhead[index]) - mean(detector.det.ylos[index])
     dz=mean(detector.det.zhead[index]) - mean(detector.det.zlos[index])
     nstep_pol=sqrt(dx^2+dy^2+dz^2)/dl
     nstep2_pol=float(round(nstep_pol*2.))
  endif
  loadct,21
  tvlct,cr,cg,cb,0,/get
  loadct,33


  for i=0,60 do tvlct,cr[i],cg[i],cb[i],i+150
  for i=0,40 do tvlct,0,100+150*i/40.,0,i+211
  for u=0,1 do begin
     if u eq 0 then begin
        index=where(detector.det.zhead lt 100.,nind) 
      ;  loadct,33
     endif else begin
        index=where(detector.det.zhead gt 100.,nind)
       ; loadct,21
     endelse
     for uu =0,nind-1 do begin
        chan=index[uu]
        ;if chan gt 2 then continue 
        if detector.det.nchan lt 100 or (chan mod 10) eq 0 then begin
           if inputs.calc_wght[0] eq 1 and inputs.ichan_wght gt 0 then begin
              if chan ne inputs.ichan_wght-1 then continue
           endif
           xhead=detector.det.xhead[chan]
           yhead=detector.det.yhead[chan]
           zhead=detector.det.zhead[chan]
           xb=xhead
           yb=yhead
           dx=xhead - detector.det.xlos[chan]
           dy=yhead - detector.det.ylos[chan]
           dz=zhead - detector.det.zlos[chan]
           if abs(zhead) gt 100. then begin ;; if a poloidal LOS
              nstep=nstep_pol
              nstep2=nstep2_pol
           endif else begin
              if zhead gt 0 then begin ;; if a toroidal LOS
                 nstep=nstep_tor
                 nstep2=nstep2_tor
              endif else begin  ;; if something else
                 ;;continue
                 nstep=sqrt(dx^2+dy^2+dz^2)/dl
                 nstep2=float(round(nstep*1.2))
              endelse
           endelse
           for ii=1.,nstep2 do begin
              dii=-ii*dx/nstep
              x=xhead+dx/dx*dii
              y=yhead+dy/dx*dii
              if zhead gt 100. then begin
                 if color_gradient then col=60./nstep2*ii+150 else col=150 
              endif else begin
                 if zhead gt 0 then begin
                     if color_gradient then col=100./nstep2*ii else col=20 
                 endif else begin
                     if color_gradient then col=211+ii/nstep2*40 else col=211 
                 endelse
              endelse
              oplot,[xb,x]/100.,[yb,y]/100.,color=col,thick=los_thick
              xb=x
              yb=y
           endfor
        endif
     endfor
  endfor
  ;; simulate the NBI source
  if 0 then begin
     isource=coords.isource
     xyz_start_src=nbgeom.xyz_src[isource,*]/100.
     for cc=0,20 do begin
        uvw_start =[ 0.d0                           $
                     ,nbgeom.bmwidra/100.* 2.d0*(randomu(seed)-0.5d0) $
                     ,nbgeom.bmwidza/100.* 2.d0*(randomu(seed)-0.5d0)]
        if uvw_start[2] lt 0. then updown=-1 else updown=1
        rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                   ,nbgeom.Brot[isource,*,*] $
                   ,nbgeom.Crot[isource,*,*],updown, xyz_start
        xyz_start=xyz_start+xyz_start_src
        ray_du=-1.d0
        focy=nbgeom.focy[isource]/100.
        focz=nbgeom.focz[isource]/100.
        ray_dv=ray_du*(uvw_start[1]/focy $
                       +tan(nbgeom.divy[0,isource]*randomn(seed)))
        ray_dw=ray_du*(uvw_start[2]/focz $
                       +tan(nbgeom.divz[0,isource]*randomn(seed)))
        uvw_ray=[ray_du,ray_dv,ray_dw]*9.     
        
        rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                   ,nbgeom.Brot[isource,*,*] $
                   ,nbgeom.Crot[isource,*,*],updown, xyz_ray
        ;;plot rays
        oplot,[xyz_start[0],xyz_start[0]+xyz_ray[0]] $
              ,[xyz_start[1],xyz_start[1]+xyz_ray[1]],thick=0.1
     endfor
     ;; plot sources
     updown=1
     uvw_start=[0.d0,1.d0,0.d0]*nbgeom.bmwidra/100.
     rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*],updown, xyz_start
     
     oplot,[ xyz_start_src[0], xyz_start_src[0]+xyz_start[0]] $
           ,[xyz_start_src[1], xyz_start_src[1]+ xyz_start[1]] $
           ,color=200
     uvw_start=[0.d0,-1.d0,0.d0]*nbgeom.bmwidra/100.
     updown=1
     rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*],updown, xyz_start
     
     oplot,[ xyz_start_src[0], xyz_start_src[0]+xyz_start[0]] $
           ,[xyz_start_src[1], xyz_start_src[1]+ xyz_start[1]] $
           ,color=200
     xyouts,[xyz_start_src[0]-3.5],[xyz_start_src[1]],'ion-source',color=200
  endif
  loadct,39
 








  ;; ----------------------------------
  ;; ----  POLOIDAL OVERVIEW ----------
  ;; ----------------------------------
  if total(plasma.denf)eq 0. then begin
     plot,[0.], /nodata, xran=rran, yran=zran $
             ,/ystyle,/xstyle,pos=pos2,/xminor,/yminor
  endif else  begin
     for i=1,254 do tvlct,256-i,256-i,256-i,i
     nr=coords.nx
     rrc=(dindgen(nr)/(nr-1.)*(max(coords.rrc_grid)-min(coords.rrc_grid)) $
          +min(coords.rrc_grid))/100.
     zplot=fltarr(nr,coords.nz)
     counter=fltarr(nr,coords.nz)*0.
     for i=0L,coords.nx-1 do begin
        for j=0,coords.ny-1 do begin
           for k=0,coords.nz-1 do begin
              dummy=min(abs(sqrt(xxc[i]^2+yyc[j]^2)-rrc),index)
              zplot[index,k]=zplot[index,k]+plasma.denf[i,j,k]*1.d6
              counter[index,k]=counter[index,k]+1.
           endfor
        endfor 
     endfor
     zplot=zplot/counter
     ;; plot colorbar
     legpos=[pos2[2],pos2[1],pos2[2]+.03,pos2[3]] <1.
     mi = 0.
     xleg = [0,1]
     expo=double(fix(alog10(ma)))
     yleg = (dindgen(255)/254.*(ma-mi)+mi)/(10.^expo)
     dumarr = fltarr(2,255)
     dumarr(0,*) = yleg
     dumarr(1,*) = yleg
     contour,dumarr,xleg,yleg $
             ,c_colors=c_colors $
             ,levels=levelvals/10.^expo $
             ,yran=[mi,ma]/10.^expo,/ysty,/fill $
             ,pos=legpos,xticks=1,yticks=1 $
             ,xtickname=(replicate(' ',2)),xticklen=-1.e-5,/xminor $
             ,ytickname=(replicate(' ',2)),yticklen=-1.e-5,/yminor 
     axis,/yaxis,yran=[mi,ma]/10.^expo,/ysty,yticks=5 $
          ,chars=chars,ytickformat='(1f3.1)'
     
     contour,zplot,rrc,zzc, xran=rran, yran=zran,/noerase $
             ,/ystyle,/xstyle,/xminor,/yminor $
             ,c_colors=c_colors $
             ,levels=levelvals,/fill,pos=pos2
   
     xyouts,pos2[2]-0.01,pos2[3]+0.015 $
            ,'fast-ion density [m!u-3!n]' $
            ,/norm,chars=chars,align=1
     xyouts,legpos[2]-.022,pos2[3]+.015,'x10!u'+strtrim(string(fix(expo)),2) $
            +'!n',/norm,chars=chars,align=0
     loadct,39
  endelse
  xyouts,pos2[0]-0.065,(pos2[1]+pos2[3])/2.,'Z [m]' $
            ,/norm,orientation=90,align=0.5
  xyouts,(pos2[0]+pos2[2])/2,0.02,'R [m]',/norm,align=0.5

  ;; PLOT DENSITY FOR PASSIVE FIDA simulation
  if keyword_set(dens) then begin
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
     counter=replicate(0.,nr,nz)
     for ix=0.,nx-1 do begin
        for iy=0.,ny-1 do begin
           dummy=min(abs(rrc-rgrid[ix,iy]),ir)
           for iz=0.,nz-1 do begin
              i=ix+coords.nx*iy+coords.nx*coords.ny*iz
              density[ir,iz]=density[ir,iz]+dens[i]
              counter[ir,iz]++
           endfor
        endfor
     endfor
     density=density/counter
     for i=1,254 do tvlct,256-i,256-i,256-i,i
     pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
     shade_surf,density,rrc,zzc,/noerase $
                ,ax=90,az=0,shades=bytscl(density,/NaN) $
                ,xtit=rtit,ytit=ztit,xran=rran,yran=zran $ 
                ,zstyle=5,/ystyle,/xstyle,pos=pos
     loadct,39
     plot,[0.], /nodata, xran=rran, yran=zran,/noerase,pos=pos $
          ,xtit='R [m]', ytit='Z [m]',/ystyle,/xstyle $
          ,xticks=5
  endif 

  ;; plot the simulation grid
  if plot_grid then oplot_fidasim_grid,coords,/pol,thick=0.5*!p.thick
  ;; plot AUG geometry
 

  ;; plot NBI
  ;if not keyword_set(dens) then begin ;; no passive fida simulation
  isource=coords.isource
  xyz_start_src=nbgeom.xyz_src[isource,*]/100.
  nii=20
  for ii=0,nii-1 do begin
     uvw_ray=[-1.d0,0.d0,0.d0]*(5.+float(ii)/float(nii-1)*4.)
     rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*] $
                ,-1,nbi_ray 
     xyz_end_src=nbgeom.xyz_src[isource,*]/100.  + nbi_ray
     r_end1=sqrt(xyz_end_src[0]^2+xyz_end_src[1]^2)
     z_end1=xyz_end_src[2]
   ;  if ii ge 1 then oplot,[r_endb1,r_end1] , [z_endb1,z_end1],thic=4,color=200
     r_endb1=r_end1
     z_endb1=z_end1
     rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*] $
                ,1,nbi_ray
     xyz_end_src=nbgeom.xyz_src[isource,*]/100.  + nbi_ray
     r_end2=sqrt(xyz_end_src[0]^2+xyz_end_src[1]^2)
     z_end2=xyz_end_src[2]

     r_mean=mean([r_end1,r_end2])
     z_mean=mean([z_end1,z_end2]) 
   if ii ge 1 then begin
        oplot,[r_meanb,r_mean] , [z_meanb,z_mean],thic=4,color=200
    ; oplot,[r_endb2,r_end2] , [z_endb2,z_end2],thic=4,color=200
  endif
     r_endb2=r_end2
     z_endb2=z_end2
     r_meanb= r_mean
     z_meanb= z_mean
  endfor
  ;; PLOT LOS
  loadct,21
  tvlct,cr,cg,cb,0,/get
  loadct,33
  for i=0,60 do tvlct,cr[i],cg[i],cb[i],i+150
  for i=0,40 do tvlct,0,100+150*i/40.,0,i+211
  for chan=0, detector.det.nchan-1 do begin
   ;if chan gt 2 then continue
     if detector.det.nchan lt 100 or (chan mod 10) eq 0 then begin
        if inputs.calc_wght[0] eq 1 and inputs.ichan_wght gt 0 then begin
           if chan ne inputs.ichan_wght-1 then continue
        endif
        xhead=detector.det.xhead[chan]
        yhead=detector.det.yhead[chan]
        zhead=detector.det.zhead[chan]
        zb=zhead
 
        rb=sqrt(xhead^2+yhead^2)
        dx=xhead - detector.det.xlos[chan]
        dy=yhead - detector.det.ylos[chan]
        dz=zhead - detector.det.zlos[chan]
        if abs(zhead) gt 100. then begin ;; if a poloidal LOS
           nstep=nstep_pol
           nstep2=nstep2_pol
        endif else begin
           if zhead gt 0 then begin ;; if a toroidal LOS
              nstep=nstep_tor
              nstep2=nstep2_tor
           endif else begin  ;; if something else
              ;continue
              nstep=sqrt(dx^2+dy^2+dz^2)/dl
              nstep2=float(round(nstep*1.2))
           endelse
        endelse
        for ii=1.,nstep2 do begin
           if zhead gt 100. then begin
              dii=-ii*dz/nstep
              x=xhead+dx/dz*dii
              y=yhead+dy/dz*dii
              z=zhead+dz/dz*dii
              if color_gradient then col=60./nstep2*ii+150 else col=150 
           endif else begin
              dii=-ii*dy/nstep
              x=xhead+dx/dy*dii
              y=yhead+dy/dy*dii
              z=zhead+dz/dy*dii
              if zhead gt 0 then begin
                 if color_gradient then col=100./nstep2*ii else col=20 
              endif else begin
                 if color_gradient then col=211+ii/nstep2*40 else col=211 
              endelse
           endelse
           r=sqrt(x^2+y^2)
           oplot,[rb,r]/100.,[zb,z]/100.,color=col,thick=los_thick
           rb=r
           zb=z
        endfor
     endif
  endfor
  if 0 then begin
     for cc=0,20 do begin
        uvw_start =[ 0.d0                           $
                     ,nbgeom.bmwidra/100.* 2.d0*(randomu(seed)-0.5d0) $
                     ,nbgeom.bmwidza/100.* 2.d0*(randomu(seed)-0.5d0)]
        if uvw_start[2] gt 0. then updown=1 else updown=-1
        rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                   ,nbgeom.Brot[isource,*,*] $
                   ,nbgeom.Crot[isource,*,*],updown, xyz_start
        xyz_start=xyz_start+xyz_start_src
        ray_du=-1.d0
        focy=nbgeom.focy[isource]/100.
        focz=nbgeom.focz[isource]/100.
        ray_dv=ray_du*(uvw_start[1]/focy $
                       +tan(nbgeom.divy[0,isource]*randomn(seed)))
        ray_dw=ray_du*(uvw_start[2]/focz $
                       +tan(nbgeom.divz[0,isource]*randomn(seed)))

        nii=400
        for ii=0,nii-1 do begin
           uvw_ray=[ray_du,ray_dv,ray_dw]*float(ii)/float(nii-1)*9.     
           
           rotate_uvw,uvw_ray,nbgeom.Arot[isource,*,*] $
                      ,nbgeom.Brot[isource,*,*] $
                      ,nbgeom.Crot[isource,*,*],updown, xyz_ray
           ;;plot rays
           r_start=sqrt(xyz_start[0]^2+xyz_start[1]^2)
           z_start=xyz_start[2]
           
           r_end=sqrt((xyz_start[0]+xyz_ray[0])^2+(xyz_start[1]+xyz_ray[1])^2)
           z_end=xyz_start[2]+xyz_ray[2]
           if ii ge 1 then oplot,[r_endb1,r_end],[z_endb1,z_end],thick=0.1
           r_endb1=r_end
           z_endb1=z_end
        endfor
     endfor
     ;; plot sources
     updown=1
     uvw_start=[0.d0,0.d0,1.d0]*nbgeom.bmwidza/100.
     rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*],updown, xyz_start

     r_start=sqrt(xyz_start_src[0]^2+xyz_start_src[1]^2)
     z_start=xyz_start_src[2]
     
     r_end=sqrt((xyz_start_src[0]+xyz_start[0])^2 $
                +(xyz_start_src[1]+ xyz_start[1])^2)
     z_end= xyz_start_src[2]+ xyz_start[2]
     oplot,[r_start,r_end],[z_start,z_end],color=200

     
     uvw_start=[0.d0,0.d0,-1.d0]*nbgeom.bmwidza/100.
     updown=-1
     rotate_uvw,uvw_start,nbgeom.Arot[isource,*,*] $
                ,nbgeom.Brot[isource,*,*] $
                ,nbgeom.Crot[isource,*,*],updown, xyz_start
     r_start=sqrt(xyz_start_src[0]^2+xyz_start_src[1]^2)
     z_start=xyz_start_src[2]
     
     r_end=sqrt((xyz_start_src[0]+xyz_start[0])^2 $
                +(xyz_start_src[1]+ xyz_start[1])^2)
     z_end= xyz_start_src[2]+ xyz_start[2]
     oplot,[r_start,r_end],[z_start,z_end],color=200
     xyouts, [r_start-1.7],[xyz_start_src[2]+0.5],'ion-source',color=200
  endif
  loadct,39
 

  if inputs.ps then begin
     file='PLOTS/inputs.eps'
     device,filename=inputs.root_dir+file
  endif


  !P.multi=0
  return 

 
end
