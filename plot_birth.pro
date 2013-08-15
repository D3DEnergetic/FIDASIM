pro plot_birth,ps=ps,path=path
  ;; PLOT ROUTINE of FIDASIM to illustrate the FAST-ION BIRTH profile of a
  ;; given NBI source
  ;; wirtten by Benedikt Geiger 2013
  if not keyword_set(path) then $
     path=dialog_pickfile(path='RESULTS/',/directory)
  runid=strsplit(path,'/',/extract,count=nid)
  runid=runid[nid-1]
  load_fidasim_results,fidasim,path
  birth_dens=fidasim.birth_dens
  ;;[fast-ions/s/grid_cell/energy_bin/pitch_bin]
  sz=size(birth_dens)
  nx=sz[1]
  ny=sz[2]
  nz=sz[3]
  nenergy=sz[4]
  npitch=sz[5]
  Efull=fidasim.inputs.einj 
  energyarr=[efull,efull/2.,efull/3.]
  pitcharr=dindgen(npitch)/(npitch)*2.-1.
  dpitch=abs(pitcharr[1]-pitcharr[0])
  nr_birth=fidasim.inputs.nr_ndmc

  if keyword_set(ps) then fidasim.inputs.ps=1
  set_plot,'X' & device, decomposed=0
  loadct,5
  !P.charthick=1. & !P.charsize=2.4 & !P.thick=1. &!P.font=-1 & !p.multi=0
  !p.background=255 & !p.color=0
  if fidasim.inputs.ps eq 1 then begin
     set_plot, 'ps'
     device,font_size=8,inches=0,/encaps $
              ,xsize=6,ysize=5,/color,bits_per_pixel=8,/helvetica
     !P.charthick=1. & !P.charsize=1. & !P.thick=1.8 &!P.font=0
     device,filename='PLOTS/birth_profile_xy_' $
            +string(fidasim.inputs.shot,f='(i5)')+'.eps'
  endif else window,0,xsize=600,ysize=500
  chars=0.9*!p.charsize
  print, 'total number of deposited fast/ions per second: ', total(birth_dens)
;;  -----------------------------------------

  contour,total(total(total(birth_dens,5),4),3),fidasim.coords.xxc $
          ,fidasim.coords.yyc,/isotropic $
          ,c_colors=(dindgen(20)+1.)*11,nlevels=20 $
          ,xtit='X [cm]',ytit='Y [cm]',/fill
  xyouts,0.3,0.8,'Birth profile',/norm
  ;; ----------------------------------------------------------
  ;; -----------  Map the BIRTH_DENS on R,Z GRID -------------
  ;; ----------------------------------------------------------
  nr=(fidasim.coords.nx+fidasim.coords.ny)/2.
  nz=fidasim.coords.nz
  rrc_grid=fidasim.coords.rrc_grid
  rrc=(dindgen(nr)/(nr-1.)*(max(rrc_grid)-min(rrc_grid))+min(rrc_grid))
  birth_dens_rz=fltarr(nr,nz,nenergy,npitch)
  for i=0L,fidasim.coords.nx-1 do begin
     for j=0,fidasim.coords.ny-1 do begin
        for k=0,fidasim.coords.nz-1 do begin
           dummy=min(abs(sqrt(fidasim.coords.xxc[i]^2 $
                              +fidasim.coords.yyc[j]^2)-rrc[*]),index)
           birth_dens_rz[index,k,*,*]+= birth_dens[i,j,k,*,*]
        endfor
     endfor 
  endfor
  ;; Conversion from [fast-ions/s/grid_cell/energy_bin/pitch_bin] to
  ;; [fast-ions/s/cm^2/energy_bin/pitch_bin]:
  dr=rrc[1]-rrc[0]
  dz=fidasim.coords.zzc[1]-fidasim.coords.zzc[0]
  birth_dens_rz=birth_dens_rz/dr/dz



  ;; --------------------------------------------------------------
  ;; ---------- Plot Birth density as a function of R an Z --------
  ;; --------------------------------------------------------------
  if fidasim.inputs.ps eq 1 then device, filename='PLOTS/birth_dens_rz_'    $
                                        + string(fidasim.inputs.shot,f='(i5)') $
                                        + '.eps' $ 
  else window,1,xsize=600,ysize=500
  ;; integrate over pitch and energy
  zplot=total(total(birth_dens_rz[*,*,*,*],4),3) ;; [fast-ions/s/cm^2]
  ma=max(zplot,/nan)
  nlevels=20
  c_colors= dindgen(nlevels)/(nlevels-1.)*253. +1.       
  levelvals=dindgen(nlevels)/(nlevels-1.)*ma
  contour,zplot,rrc,fidasim.coords.zzc $
          ,c_colors=c_colors,levels=levelvals,/fill  $
          ,xran=[100,250],yran=[-100,100],/isotropic $
          ,xtit='R [cm]',ytit='Z [cm]'
  ;; PLOT LEGEND
  legpos=[!x.window[1],!y.window[0],!x.window[1]+0.04,!y.window[1]]
  xleg = [0,1]
  expo=double(fix(alog10(ma)))
  yleg = (dindgen(255)/254.*ma)/(10.^expo)
  dumarr = fltarr(2,255)
  dumarr(0,*) = yleg
  dumarr(1,*) = yleg
  contour,dumarr,xleg,yleg,/noerase $
          ,c_colors=c_colors,/fill $
          ,levels=levelvals/10.^expo $
          ,yran=[0.,ma]/10.^expo,/ysty $
          ,pos=legpos,xticks=1,yticks=1 $
          ,ytickname=(replicate(' ',2)),yticklen=-1.e-5,yminor=1 $
          ,xtickname=(replicate(' ',2)),xticklen=-1.e-5,xminor=1
  axis,/yaxis,yran=[0.,ma]/10.^expo,/ysty,yticks=5 $
       ,chars=chars,ytickformat='(1f3.1)'
  xyouts,legpos[2],legpos[3]+0.01,'deposited fast-ions [10!u' $
         +strtrim(string(fix(expo)),2)+'!n/s/cm!u2!n]' $
         ,chars=chars,/norm,/align
 


  ;; --------------------------------------------------------------
  ;; ------------ Plot Birth density as a function of R  ----------
  ;; --------------------------------------------------------------
  if fidasim.inputs.ps eq 1 then device, filename='PLOTS/birth_dens_r_'    $
                                         +string(fidasim.inputs.shot,f='(i5)') $
                                         +'.eps' $
  else window,2,xsize=600,ysize=500
  zplot=total(zplot,2)*dz ;; integrate along z
  ma=max(zplot)
  expo=double(fix(alog10(ma)))
  plot,[0.],/nodata,xran=[100,250],yran=[0,ma/10^expo] $
       ,xtit='R [cm]',ytit='10!u'+strtrim(string(fix(expo)),2) $
       +'!nfast-ions/(s cm)' 
  oplot,rrc,zplot/10^expo


  ;; ----------------------------------------------------------
  ;; ------  CALCULATE RANDOM BIRTH POSITIONS of markers ------
  ;; ----------------------------------------------------------
  for i=0,2000 do begin
     birth_dens2=birth_dens/total(birth_dens) * nr_birth*(1.+i*0.01)
     if total(birth_dens2) gt nr_birth then goto, here
  endfor
  stop
  here:
  tb=total(birth_dens2)+100
  birth_pos=fltarr(3,tb)
  birth_p=fltarr(tb)
  birth_e=fltarr(tb)
  cc=0L
  for i=0,nx-1 do begin
     for j=0,ny-1 do begin   
        for k=0,nz-1 do begin
           for ie=0,nenergy-1 do begin
              zz=reform(birth_dens2[i,j,k,ie,*])
              if total(zz) gt 0 then begin
                 nib=total(zz)
                 ben=nib mod 1
                 nib=floor(nib)
                 if ben gt randomu(seed) then nib=floor(nib)+1.

                 for ii=0L,nib-1 do begin
                    birth_pos[*,cc]=[fidasim.coords.xx[i]   $
                                     ,fidasim.coords.yy[j]  $
                                     ,fidasim.coords.zz[k]] $
                                    +[fidasim.coords.dx $
                                      ,fidasim.coords.dy $
                                      ,fidasim.coords.dz] $
                                    *randomu(seed,3)
                    if nib eq 1 then begin
                       dummy=max(zz,ipitch)
                    endif else begin
                       ww=zz/max(zz)
                       for pp=0,1000 do begin
                          ran=randomu(seed,2)
                          ipitch=round(ran[0]*(npitch-1.))
                          if ww[ipitch] gt ran[1] then goto, out
                       endfor
                       print,'problem with pitch'
                       dummy=max(zz,ipitch)
                       out:
                    endelse
                    birth_p[cc]=pitcharr[ipitch]+dpitch*randomu(seed,1)
                    birth_e[cc]=energyarr[ie]
                    cc++
                 endfor ;; loop over nib

              endif
           endfor
        endfor
     endfor
  endfor
  birth_pos=birth_pos[*,0:cc-1]
  birth_e=birth_e[0:cc-1]
  birth_p=birth_p[0:cc-1]
  print, cc ,' fast-ion vectors have been calculated!!'


  ;; convert to RZ-PHI!!
  birth_rzphi=fltarr(3,cc)
  birth_rzphi[0,*]=sqrt(birth_pos[0,*]^2+birth_pos[1,*]^2)
  birth_rzphi[1,*]=birth_pos[2,*]
  birth_rzphi[2,*]=atan(birth_pos[1,*]/birth_pos[0,*]) 
  ;; CORRECT THE ANGLE FOR NEGATIVE X and Y values!
  ;; if x<0 and y>0: +pi
  index=where(birth_pos[0,*] lt 0 and birth_pos[1,*] gt 0,nind)
  if nind ne 0 then  birth_rzphi[2,index]=birth_rzphi[2,index]+!pi
  ;; if x<0 and y<0: +pi
  index=where(birth_pos[0,*] lt 0 and birth_pos[1,*] lt 0,nind)
  if nind ne 0 then birth_rzphi[2,index]=birth_rzphi[2,index]+!pi
  ;; if x>0 and y<0: +2pi
  index=where(birth_pos[0,*] gt 0 and birth_pos[1,*] lt 0,nind)
  if nind ne 0 then birth_rzphi[2,index]=birth_rzphi[2,index]+2.*!pi



  ;; --------------------------------------------------
  ;; --------- PLOT THE RANDOM BIRTH PROFILE!! --------
  ;; --------------------------------------------------
  ;; NOW PLOT THE RZ-BIRTH PROFILE
  if fidasim.inputs.ps eq 1 then device, filename='PLOTS/birth_profile_rz_'    $
                                        + string(fidasim.inputs.shot,f='(i5)') $
                                        + '.eps' $ 
  else window,3,xsize=600,ysize=500
  
  plot,[0.],/nodata,xran=[100,250],yran=[-100,100],/isotropic $
       ,xtit='R [cm]',ytit='Z [cm]'
  index=where(birth_p lt -.5,nind)
  if nind gt 0 then $
     oplot, [birth_rzphi[0,index]],[birth_rzphi[1,index]],psym=2 $
            ,col=254,symsize=0.5
  xyouts,0.25,0.3+0.05,'pitch < -0.5',color=254,/norm
  index=where(birth_p lt -0.4 and birth_p gt -0.5,nind)
   if nind gt 0 then $
      oplot, [birth_rzphi[0,index]],[birth_rzphi[1,index]],psym=4 $
             ,color=20,symsize=0.5
   xyouts,0.25,0.3,'-0.5 < pitch < -0.4',color=20,/norm
   index=where(birth_p gt -0.4,nind)
   if nind gt 0 then $
      oplot, [birth_rzphi[0,index]],[birth_rzphi[1,index]],psym=5 $
             ,color=60,symsize=0.5
   xyouts,0.25,0.3-0.05,'-0.4 < pitch',color=60,/norm


  ;; -------------------------------------------------------
  ;; --------- WRITE THE BIRTH PROFILE to a file !! --------
  ;; -------------------------------------------------------
  file = path+'/birth_profile.dat'
  openw, 55, file
  printf,55,'# FIDASIM birth profile created: ', systime(),' Version 1.0'
  printf,55, fidasim.inputs.shot         ,f='(i6,"         # shotnumber")'  
  printf,55, fidasim.inputs.time,f='(1f8.5,"       # time")'
  printf,55, cc,f='(i8,"       # fast-ion markers")'
  printf,55,'#     R [cm],     Z [cm],      phi [rad],   energy [keV],   pitch'
  for i=0L,cc-1 do begin
     printf,55, birth_rzphi[0,i], birth_rzphi[1,i], birth_rzphi[2,i] $
            , birth_e[i], birth_p[i]
  endfor
  close,55
  if fidasim.inputs.ps eq 1 then device,/close
  stop
end                             ;of programm

