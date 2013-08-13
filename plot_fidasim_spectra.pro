pro plot2d,zz,xx,yy,emax,xmar,thick,transp=transp,zmax,pos
  shade_surf,zz,xx,yy,$ 
             ax=90,az=0,shades=bytscl(zz,top=253,/NaN,max=zmax),$
             xran=[0.,emax],yran=[1.,-1.], zstyle=5, /ystyle,$
             xtit='Energy [keV]',ytit='Pitch',xmar=xmar
  legpos=[0.95,!y.window[0],0.99,!y.window[1]]
  pos=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
  mi = min(zz,/NaN)
  ma = zmax
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
     ytickname(i) = string(yvalue,f='(1e7.1)')      
  endfor
  tvlct,255,255,255,0
  shade_surf,dumarr,xleg,yleg,/noerase, $
             ax=90,az=0,shades=bytscl(dumarr,top=253,/NaN),$
             yran=[mi,ma],zstyle=5,/ystyle,$
             position=legpos,xticks=1,yticks=legticks-1,chars=.7,$
             xticklen=0.0001,yticklen=-0.01,$
             xtickname=(replicate(' ',2)),$
             ytickname=ytickname,xthick=thick,ythick=thick
  tvlct,0,0,0,0
  plot,[0.],/nodata $
       ,xran=[0.,emax],yran=[1.,-1.], zstyle=5, /ystyle $
       ,ytit='Pitch',/noerase,xthick=thick $
       ,ythick=thick,pos=pos
end

pro oplot_transp,zz,xx,yy,emax,thick,pos
  tvlct,155,155,155,0
  contour,zz,xx,yy,/noerase, $ 
          xran=[0.,emax],yran=[1.,-1.], zstyle=5, /ystyle,$
          color=0,pos=pos, nlevels=8
  tvlct,255,255,255,0
  plot,[0.],/nodata $
       ,xran=[0.,emax],yran=[1.,-1.], zstyle=5, /ystyle,$
       /noerase,xthick=thick $
       , ythick=thick,pos=pos
end


pro plot_fidasim_spectra,ps=ps,product=product,runid=runid
 ;; PLOT SETTINGS
  linthick=1.
  loadct,39    ;green=150;blue=50;yellow=200;red=254;black=0;white=255
  set_plot,'X' & device, decomposed=0
  !P.charthick=1. & !P.charsize=1.5 &!p.font=-1 & !P.thick=1.0
  if keyword_set(ps) then begin 
     set_plot, 'ps'
     directory='PLOTS/'
     device, font_size=8, inches=0 , /encaps $
             ,xsize=7.9, ysize=6, /color, bits_per_pixel=8 
     device, /helvetica,font_index=3
     device, /symbol,font_index=4                                    
     device, /helvetica
     linthick=1.
     linthick=1
     !P.charsize=1. & !P.charthick=3. & !P.thick=1.5 & !P.font=0
  endif    
  !P.background=255 & !P.color=254 
  
  
  
  ;; LOAD FIDASIM DATA 
  if keyword_set(runid) then path=runid else $
     path=dialog_pickfile(path='RESULTS/',/directory)
  load_fidasim_results, fidasim,path
  shot=fidasim.inputs.shot
  time=fidasim.inputs.time


  if fidasim.inputs.calc_wght eq 1 then begin
     wfunct=fidasim.wfunct.weight_tot
     cdf_file=fidasim.inputs.transp_runid
     print, 'cdf file: ' ,cdf_file
     print, 'nr weight:', fidasim.wfunct.nen
     
     
     
     
     ;; LOAD transp fast-ion distribution function
     cdfid=NCDF_Open(cdf_file,/nowrite)
     ;;Retrieve signals
     ;;--------------------------------------------------------
     ncdf_varget, cdfid,'TRANSP_RUNID', runid
     ncdf_varget, cdfid,'TIME' , cdf_time    
     print, runid,'cdf time:', cdf_time
     ncdf_varget, cdfid,'R2D'  , r2d             ; rposition of cells
     ncdf_varget, cdfid,'Z2D'  , z2d             ;zposition of cells
     ncdf_varget, cdfid,'E_D_NBI', energyarr_transp ; central values
     ncdf_varget, cdfid,'A_D_NBI', pitcharr_transp  ; central values
     ncdf_varget, cdfid,'NTOT_D_NBI',ntot           ; total number of fast ions 
     nzones_id=ncdf_varid(cdfid,'N_ZONEROWS')
     ncdf_varget,cdfid,nzones_id,n_zones
     print,'NZONE_FB:',n_zones
     ncdf_varget, cdfid,'F_D_NBI', FBM ; beam distribution function
     ;; Fetch selected species
     ;; CLF: need "spc_lab" to generate file name
     spc=1
     spc_id=ncdf_varid(cdfid,'SPECIES_'+strtrim(string(spc),1))
     ncdf_varget,cdfid,spc_id,spc_lab
     print, 'spc_lab:', spc_lab
     spc_lab=strtrim(string(spc_lab))
     IF (N_ELEMENTS(fus) GT 0) THEN spc_lab=spc_lab+'_FUSN'
     ;; Fetch pitch angle grid and grid boundary values
     a_d_id=ncdf_varid(cdfid,'A_'+spc_lab)
     ncdf_varget,cdfid,a_d_id,a_d
     ;; n_sym=1 (updown symmetric equilibrium), 2 (asymmetric)
     symflag_id=ncdf_varid(cdfid,'SYMFLAG')
     ncdf_varget,cdfid,symflag_id,symflag
     if (symflag EQ 0) THEN BEGIN 
        n_sym=2
     ENDIF ELSE BEGIN
        n_sym=1
     ENDELSE
     NCDF_Close,cdfid
     ;;====================
     ;;passing vs. trapped
     ;;====================
     n_E=(SIZE(fbm))(1)
     n_pit=(SIZE(fbm))(2)
     n_cells=(SIZE(fbm))(3)
     ;; same rho_lab for same rho value, takes values 0...n_zones-1
     rho_lab=intarr(n_cells)
     FOR j_rho=1,n_zones DO BEGIN
        min_rho=n_sym*j_rho*(j_rho-1)
        max_rho=n_sym*j_rho*(j_rho+1)
        rho_lab(min_rho:max_rho-1)=j_rho-1
     ENDFOR
     ;; Amount of theta points at a given rho
     the_size=fltarr(n_zones)
     FOR j_rho=0,n_zones-1 DO BEGIN
        the_size(j_rho)=2*n_sym*(j_rho+1)
     END
     rmaj_min=fltarr(n_cells)
     count=0
     FOR j_rho=0,n_zones-1 DO BEGIN
        rmaj_min[j_rho]=MIN(r2d[count:count+the_size(j_rho)-1])
        count=count+the_size(j_rho)
     END
     trap_pit=fltarr(n_cells)
     trap_pit[*]=1.-rmaj_min[rho_lab(*)]/r2d[*]
     ;;----------- Check the time
     print, ' CDF file time:',cdf_time  
     ;;-------------Convert eV/solid angel/4pi-> keV/dpitch
     energyarr_transp=energyarr_transp*1.0d-3  ;; fidasim needs energy in kev  
     fbm=fbm*1.0d3*0.5   
     nen_transp=n_elements(energyarr_transp)
     npitch_transp=n_elements(pitcharr_transp) 
     ;;reverse array
     pos=npitch_transp-1-indgen(npitch_transp)
     fbm=reform(fbm[*,pos,*])  
     ;; define dE and dptich
     dE_transp=energyarr_transp[1]-energyarr_transp[0]
     dpitch_transp=abs(pitcharr_transp[1]-pitcharr_transp[0])

     radiance= fltarr(fidasim.wfunct.nchan,fidasim.wfunct.nwav)
     fbm_tot = fltarr(fidasim.wfunct.nchan,fidasim.wfunct.nwav)
     mean_fbm= fltarr(nen_transp,npitch_transp,fidasim.wfunct.nchan)
  endif
  
  
  
  
  
  
  ;; ------------------------------------------
  ;; START LOOP over the different channels----
  ;; ------------------------------------------
  for ichan=0, fidasim.los.nchan-1 do begin
     print, 'ichan:', ichan
     rlos=sqrt(fidasim.los.xyzlos[ichan,0]^2  + $
               fidasim.los.xyzlos[ichan,1]^2)/100. ;m
     zlos=fidasim.los.xyzlos[ichan,2]/100.         ;[m]
     print, 'R= ',strtrim(rlos,2), ' m'
     print, 'Z= ',strtrim(zlos,2), ' m'
     
     if fidasim.los.xyzhead[ichan,2] lt 100. then kind='toroidal'else $
        kind='poloidal'
     
     
     if fidasim.inputs.calc_wght eq 1 then begin
        if fidasim.wfunct.ichan ne ichan and fidasim.wfunct.ichan gt 0 then continue

        ;; ------------------------------------------
        ;;------------ CALCUATE mean value of fbm----
        ;; ------------------------------------------
        trapped_pitch=0.d0
        rad=0.d0
        for i=0,fidasim.coords.nx-1 do begin
           for j=0,fidasim.coords.ny-1 do begin
              for k=0,fidasim.coords.nz-1 do begin
                 los_wght=fidasim.los.los_weight[i,j,k,ichan]
                 if los_wght gt 0. then begin
                    ;; determine mean values like the halo density along LOS
                    wght=(  fidasim.neutrals.fdens[i,j,k,2] + $
                            fidasim.neutrals.hdens[i,j,k,2] + $
                            fidasim.neutrals.tdens[i,j,k,2] + $
                            fidasim.neutrals.halodens[i,j,k,2]) * los_wght
                    rad=rad+wght
                    
                    rrc=sqrt(fidasim.coords.xxc[i]^2+fidasim.coords.yyc[j]^2)
                    zzc=fidasim.coords.zzc[k]
                    dr=2.       ;[cm]
                    dz=2.       ;[cm]
                    dummy=min((r2d-rrc+dr)^2+(z2d-zzc+dz)^2,fbm_index1)
                    dummy=min((r2d-rrc-dr)^2+(z2d-zzc+dz)^2,fbm_index2)
                    dummy=min((r2d-rrc+dr)^2+(z2d-zzc-dz)^2,fbm_index3)
                    dummy=min((r2d-rrc-dr)^2+(z2d-zzc-dz)^2,fbm_index4)
                    dummy=min((r2d-rrc)^2+(z2d-zzc)^2,fbm_index5)
                    mean_fbm[*,*,ichan] = mean_fbm[*,*,ichan] + ( $
                                       fbm[*,*,fbm_index1]+fbm[*,*,fbm_index2] + $
                                       fbm[*,*,fbm_index3]+fbm[*,*,fbm_index4] + $
                                       fbm[*,*,fbm_index5])/5.*wght  
                    trapped_pitch = trapped_pitch +( $
                                    trap_pit[fbm_index1]+trap_pit[fbm_index2] + $
                                    trap_pit[fbm_index3]+trap_pit[fbm_index4] + $
                                    trap_pit[fbm_index5])/5.*wght
                 endif
              endfor
           endfor
        endfor
        mean_fbm[*,*,ichan]=mean_fbm[*,*,ichan]/rad
        trapped_pitch=trapped_pitch/rad
        print, 'trapped_pitch :', trapped_pitch
        
        
        ;;------------------------------------------------------------
        ;; map FBM on the energy and pitch grid of the weight function
        ;;------------------------------------------------------------
        mean_fbm2=fltarr(fidasim.wfunct.nen,fidasim.wfunct.npitch)
        for ie=0,fidasim.wfunct.nen-1 do begin
           dummy=min(abs(energyarr_transp-fidasim.wfunct.energyarr[ie]),eindex)
           for ip=0,fidasim.wfunct.npitch-1 do begin
              dummy=min(abs(pitcharr_transp-fidasim.wfunct.pitcharr[ip]),pindex)
              mean_fbm2[ie,ip]=mean_fbm[eindex,pindex,ichan]
              ;;[fast_ion/cm^3/dP/dE]
           endfor
        endfor
        
        
        ;;------------------------------------------------------------
        ;; ------ plot the weight functions -------------------------
        ;;------------------------------------------------------------
        if not keyword_set(output) then begin
           ;; plot weight functions
           loadct,39
           for i=0,60 do tvlct,255-i*4.,255-i*4.,255,i
           tvlct,0,0,0,254
           xmar=[7,9]
           xx=fidasim.wfunct.energyarr
           yy=fidasim.wfunct.pitcharr
           dwav=fidasim.wfunct.dwav
           emax=fidasim.wfunct.emax
           
           
           ;; if not keyword_set(ps) then window,1
           print, 'press enter to continue, press "S" to skip the plots'
           for ii=0,fidasim.wfunct.nwav-1 do begin
              ;;  print, ii
              ;;  if ii ne 100 then continue
              if not keyword_set(product) then begin
                 if keyword_set(ps) then device,filename='PLOTS/'+kind $
                                                +'_wfunct_chan' $
                                                +strtrim(string(ichan+1,f='(i2)'),2)  $
                                                +'_'+strtrim(string(ii+1),2) $
                                                +'.eps' else  wset,0
                 
                 ;; plot weight function
                 zz= wfunct[ichan,ii,*,*]
                 zmax=max(zz)
                 tvlct,255,255,255,0
                 plot2d,zz,xx,yy,emax,xmar,thick,zmax,pos
                 oplot_transp,mean_fbm[*,*,ichan],energyarr_transp $
                              ,pitcharr_transp,emax,thick,pos
                 xyouts,0.13,0.94,'#'+strtrim(string(shot),2) $
                        +'@'+string(time,f='(1f5.3)') $
                        +'s, R='+string(rlos,f='(1f6.2)')+'m',/norm,chars=.8
                 xyouts,0.99,0.94 , '[ph*cm/(s*fast_ion)]' $
                        ,/norm,chars=.8,alignment=1.
                 xyouts,0.35,0.25,'wavelength [nm]: ' $
                        + string(fidasim.wfunct.central_wavel[ii] $
                                 -0.5*dwav,f='(1f5.1)') $
                        +' - '+string(fidasim.wfunct.central_wavel[ii]+0.5*dwav $
                                      ,f='(1f5.1)'),/norm,chars=.9
                 
                 ;; plot trapped-passing boundary
                 tvlct,155,155,155,0
                 oplot, [0,300],[trapped_pitch,trapped_pitch],linesty=1,col=0
                 oplot, [0,300],[-trapped_pitch,-trapped_pitch],linesty=1,col=0
                 tvlct,255,255,255,0
                 wait, 0.01
                 dummy='' 
                 
              endif else begin
                 ;; plot the product of the weight function and FBM_mean2
                 if keyword_set(ps) then begin
                    device,filename='PLOTS/'+kind+'_prod_chan' $
                           +strtrim(string(ichan+1,f='(i2)'),2) $
                           +'_'+strtrim(string(ii+1),2) $
                           +'.eps'
                 endif else  wset,0
                 prod=replicate(0.,fidasim.wfunct.nen,fidasim.wfunct.npitch)
                 for ie=0,fidasim.wfunct.nen-1 do begin
                    for ip=0,fidasim.wfunct.npitch-1 do begin
                       prod[ie,ip]=mean_fbm2[ie,ip]*wfunct[ichan,ii,ie,ip]
                       ;;--> [ph/(s cm^2 dP keV)]
                    endfor
                 endfor
                 zz=prod
                 zmax=max(zz)
                 tvlct,255,255,255,0
                 plot2d,zz,xx,yy,emax,xmar,thick,zmax,pos
                 oplot_transp,mean_fbm[*,*,ichan],energyarr_transp $
                              ,pitcharr_transp,emax,thick,pos
                 xyouts,0.13,0.94,'#'+strtrim(string(shot),2) $
                        +'@'+string(time,f='(1f5.3)') $
                        +'s, R='+string(rlos,f='(1f6.2)')+'m',/norm,chars=.8
                 xyouts,0.99,0.94 , '[ph*cm/(s*fast_ion)]' $
                        ,/norm,chars=.8,alignment=1.
                 xyouts,0.35,0.25,'wavelength [nm]: ' $
                        + string(fidasim.wfunct.central_wavel[ii] $
                                 -0.5*dwav,f='(1f5.1)') $
                        +' - '+string(fidasim.wfunct.central_wavel[ii]+0.5*dwav $
                                      ,f='(1f5.1)'),/norm,chars=.9
                 
                 ;; plot trapped-passing boundary
                 tvlct,155,155,155,0
                 oplot, [0,300],[trapped_pitch,trapped_pitch],linesty=1,col=0
                 oplot, [0,300],[-trapped_pitch,-trapped_pitch],linesty=1,col=0
                 tvlct,255,255,255,0
                 wait, 0.01
                 dummy='' 
              endelse ;; plot product
              if not keyword_set(ps) then read,dummy
              if dummy eq 's' or dummy eq 'S' then goto, noplot
           endfor ;; loop over wavelength
        endif     ;;plot the weight function
        noplot:
        ;;------------------------------------------------------------
        ;; ------ CALCULATE SYNTHETIC SPECTRA AND PROFIELES ----------
        ;;------------------------------------------------------------
        for ii=0,fidasim.wfunct.nwav-1 do begin
           ;; PRODUCT with fast-ion distribution funciton
           prod=replicate(0.,fidasim.wfunct.nen,fidasim.wfunct.npitch)
           for ie=0,fidasim.wfunct.nen-1 do begin
              for ip=0,fidasim.wfunct.npitch-1 do begin
                 prod[ie,ip]=mean_fbm2[ie,ip]*wfunct[ichan,ii,ie,ip]
                 ;;--> [ph/(s cm^2 dP keV)]
              endfor
           endfor
           radiance[ichan,ii]=total(prod[*,*])*fidasim.wfunct.dE    *  $
                           fidasim.wfunct.dpitch/(4.d0*!pi)*1.d4 /  $
                           fidasim.wfunct.dwav
           ;;--> [ph/(s m^2 sr nm)]
        endfor
     endif

     loadct,39 ;; plot spectrum  
    
     if fidasim.los.xyzhead[ichan,2] lt 100. then kind='toroidal'else $
        kind='poloidal'
     if keyword_set(ps) then  $
        device,filename='PLOTS/' $
               +kind+'_spectrum_weight_function_' $
               +'chan'+strtrim(string(ichan+1,f='(i2)'),2) $
               +'.eps' else window,3
     xran=[647,665]
     yran=[[1.e15],1.e18]
     plot, [0.],/nodata,xran=xran,yran=yran,/xsty $
           , ytit='Intensity [Ph/(s m^2 nm sr)]', xtit='lambda [nm]' $
           ,xthick=linthick,ythick=linthick,/ylog,color=0

     brems=fidasim.spec.brems[*,ichan]
     oplot,fidasim.spec.lambda,fidasim.spec.fida[*,ichan]+brems ,color=254
     oplot,fidasim.spec.lambda,fidasim.spec.halo[*,ichan]+brems ,color=150
     oplot,fidasim.spec.lambda,fidasim.spec.full[*,ichan]+brems,color=210
     oplot,fidasim.spec.lambda,fidasim.spec.half[*,ichan]+brems ,color=200
     oplot,fidasim.spec.lambda,fidasim.spec.third[*,ichan]+brems ,color=190
 

     if fidasim.inputs.calc_wght eq 1 then begin
        if fidasim.wfunct.ichan gt 0 and ichan eq 0 then begin
           ichan=fidasim.wfunct.ichan
        endif else ichan=ichan
        oplot,fidasim.wfunct.central_wavel,radiance[ichan,*]+brems ,color=50
        xyouts,0.22,0.7,'-  weight funciton method', /norm, color=50,align=0
     endif
     xyouts,0.22,0.8,'Synthetic FIDA spectrum by:', /norm, color=0,align=0
     xyouts,0.22,0.75,'-  Monte Carlo approach', /norm, color=254,align=0
     
     rlos=sqrt(fidasim.los.xyzlos[ichan,0]^2  + $
               fidasim.los.xyzlos[ichan,1]^2)/100. ;m
     xyouts,0.2,0.87,'#'+strtrim(string(shot),2) $
            +'@'+string(time,f='(1f5.3)') $
            +'s, R='+string(rlos,f='(1f6.2)')+'m',/norm,chars=.8,color=0
     xyouts,0.42,0.2,kind+' LOS',/norm,color=0
     
     dummy=''
     if not keyword_set(ps) then read,dummy
     
  endfor   ;; loop over channels
  wspec=fidasim.wfunct.central_wavel
  spec=fidasim.spec
  save,spec,radiance,wspec,filename='spec.sav'
  
  if keyword_set(ps) then device, /close
end






