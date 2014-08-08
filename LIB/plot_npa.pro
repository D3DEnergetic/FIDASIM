pro plot_npa,histo,energy_arr,pitch_arr,distri,ps=ps,path=path,chan=chan,currentmode=currentmode
  ;; ROUTINE OF FIDASIM to plot resulting NPA fluxes
  ;; written by Philipp Scheider and Benedikt Geiger 2013
  ;; plot settings

  set_plot,'X' & device, decomposed=0
  !P.charthick=1. & !P.charsize=1.5 & !P.thick=1.0 
  loadct,39
  !P.background=255 & !P.color=0 
  if keyword_set(ps) then begin  
     set_plot, 'ps',/COPY
     device, filename='PLOTS/npa_top_down.ps',font_size=8, inches=0 , encaps=1 $
             ,xsize=5, ysize=4.5, /color, bits_per_pixel=8, /helvetica
     !P.charsize=1. & !P.charthick=1. & !P.thick=3.5 & !P.font=0
  endif

  ;; load fidasim results
  if not KEYWORD_SET(path) then begin
     path=DIALOG_PICKFILE(path='RESULTS/',/directory)
     runid=strsplit(path,'/',/extract,count=nid)
     runid=runid[nid-1]
  endif
  if not keyword_set(chan) then chan=0
  load_results, path,fidasim
  shot=fidasim.inputs.shot
  time=fidasim.inputs.time
  los=fidasim.los
  ;; density along NBI path
  dens=total(fidasim.neutrals.halodens[*,*,*,*],4) $
       +total(fidasim.neutrals.fdens[*,*,*,*],4)   $
       +total(fidasim.neutrals.hdens[*,*,*,*],4)   $
       +total(fidasim.neutrals.tdens[*,*,*,*],4)   

  ;; initial position where fast-neutral is born (ion is neutraliszed)
  npaipos=fidasim.npa.ipos[*,*,chan]

  ;; position where neutral is detected!
  npafpos=fidasim.npa.fpos[*,*,chan]

  ;; velocity vector of fast-neutral
  npav=fidasim.npa.v[*,*,chan]

  ;; weight (particle number per marker)
  npawght=fidasim.npa.wght[*,chan]

  flux=fidasim.npa.flux[*,chan]
  energy_arr=fidasim.npa.energy

  ;; npa weight: Particles/s/cm^2
  ww=where(npawght ne 0,nnpa)

  npawght=npawght[ww]
  npav=npav[ww,*]
  npaipos=npaipos[ww,*]
  npafpos=npafpos[ww,*]

  if keyword_set(currentmode) then begin
    ; Use simple model based on Fig. 2 of Shinohara et al., RSI 75 (2004) 3640
    eintercept=25.        ; keV--value depends on foil & detector
    cmode=0.
  end

  ;; plot initial and end position of fast-ion trajectories on
  ;; top-down view
  contour,total(dens,3),fidasim.grid.x_grid[*,*,0],fidasim.grid.y_grid[*,*,0] $
          ,c_colors=(dindgen(20)+1.)*11,nlevels=20 $
          ,/isotropic,xtit='X [cm]',ytit='Y [cm]'
  oplot,fidasim.grid.x_grid[*,*,0],fidasim.grid.y_grid[*,*,0],psym=3
  for ichan=0,los.nchan-1 do begin
      if fidasim.los.chan_id[ichan] ne 1 then continue
      oplot,[los.xlens[ichan],los.xlos[ichan]] $
          ,[los.ylens[ichan],los.ylos[ichan]]
  endfor
  for i=0,nnpa-1 do oplot,[npaipos[i,0]],[npaipos[i,1]],psym=3,color=254./nnpa*i

  ;; PLOT Histogram
  if keyword_set(ps) then device, filename='PLOTS/npa_histogram.eps' $
  else window,1

  ;; Define arrays for histograms
  ;; energy array:
  nen  = n_elements(energy_arr)
  dE   = energy_arr[1]-energy_arr[0]

  ;;pitch array:
  npitch = 40
  pmin   = -1.
  pmax   = 1.
  pitch_arr = pmin+dindgen(npitch)/(npitch-1.)*(pmax-pmin)
  dpitch = abs(pitch_arr[1]-pitch_arr[0])

  ;; rho_array
  nrhop = 56
  rmin  = 0.0
  rmax  = 1.1
  rhop_arr = rmin+DINDGEN(nrhop)/(nrhop-1)*(rmax-rmin)

  ;; detector heads (here only 1)
  ndet=1
  distri=fltarr(ndet,nen,npitch)

; loop over all particles which reach the detector
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
  mass=2.0*1.66e-27                                 ; [kg]
  ec=1.6e-19   
  !P.multi = [0,0,2]
  pEnergy = DBLARR(nnpa)

  for idet=0,ndet-1 do begin ;; loop over detectors
     for i=0L,nnpa-1 do begin  ;; loop over markers that reach the detector
        ;; velocity
        vn=npav[i,*]
        vabs=sqrt(vn[0]^2+vn[1]^2+vn[2]^2)
        vec=vn/vabs

        ;; calculate energy of neutral
        energy=0.5*mass*vabs^2/(ec*1.e3)*1.e-4 ; [keV]

        if keyword_set(currentmode) then begin
          if energy gt eintercept then cmode+=(energy-eintercept)*npawght[i]
        endif
        ;; get magnetic field where neutral was born
        dummy=min(abs(fidasim.grid.x_grid[*,0,0]-npaipos[i,0]),xind)
        dummy=min(abs(fidasim.grid.y_grid[0,*,0]-npaipos[i,1]),yind)
        dummy=min(abs(fidasim.grid.z_grid[0,0,*]-npaipos[i,2]),zind)

		bvec=[fidasim.plasma.bx[xind,yind,zind],$
		      fidasim.plasma.by[xind,yind,zind],$
			  fidasim.plasma.bz[xind,yind,zind]]

        ;;normalize
        bvec=bvec/sqrt(bvec[0]^2+bvec[1]^2+bvec[2]^2)

        ;; calculate pitch of neutral
        pitch=(vec[0]*bvec[0]+vec[1]*bvec[1]+vec[2]*bvec[2]) $
              *fidasim.plasma.btipsign
        dummy=min(abs(energy_arr-energy),eindex)
        dummy=min(abs(pitch_arr-pitch),pindex)    

        if keyword_set(currentmode) then begin
            cm=0
            if energy gt eintercept then cm=(energy-eintercept)
	        distri[idet,eindex,pindex]=distri[idet,eindex,pindex] + cm
		endif else begin    
    		distri[idet,eindex,pindex]=distri[idet,eindex,pindex]+npawght[i]
        endelse
        pEnergy[i] = eindex
     endfor

     contour,distri[idet,*,*],energy_arr,pitch_arr $
             ,c_colors=indgen(20)*12,nlevels=20,/fill,yran=[-1,1] $
             ,ytit='Pitch',xtit='Energy [keV]'
     plot, energy_arr, flux,ytit='neutrals/s',xtit='Energy [keV]',psym=10

  endfor
  if keyword_set(ps) then device, /close

  if keyword_set(currentmode) then print,'Current mode:',cmode

  !p.multi=0
end 

