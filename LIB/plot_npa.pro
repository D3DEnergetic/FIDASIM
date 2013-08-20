pro plot_npa,ps=ps,path=path
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

  load_results, path,fidasim
  shot=fidasim.inputs.shot
  time=fidasim.inputs.time

  ;; density along NBI path
  dens=total(fidasim.neutrals.halodens[*,*,*,*],4) $
       +total(fidasim.neutrals.fdens[*,*,*,*],4)   $
       +total(fidasim.neutrals.hdens[*,*,*,*],4)   $
       +total(fidasim.neutrals.tdens[*,*,*,*],4)   

  ;; initial position where fast-neutral is born (ion is neutraliszed)
  npaipos=fidasim.npa.npaipos

  ;; position where neutral is detected!
  npafpos=fidasim.npa.npafpos

  ;; velocity vector of fast-neutral
  npav=fidasim.npa.npav

  ;; weight (particle number per marker)
  npawght=fidasim.npa.npawght

  ;; npa weight: Particles/s/cm^2
  headsize=mean(los.headsize)
  npawght=npawght*(!pi*headsize^2)
  nnpa=n_elements(npawght)

  ;; plot initial and end position of fast-ion trajectories on
  ;; top-down view
  contour,total(dens,3),fidasim.grid.x_grid[*,*,0],fidasim.grid.y_grid[*,*,0] $
          ,c_colors=(dindgen(20)+1.)*11,nlevels=20 $
          ,/isotropic,xtit='X [cm]',ytit='Y [cm]'

  for ichan=0,los.nchan-1 do oplot,[los.xlens[ichan],los.xlos[ichan]] $
                                   ,[los.ylens[ichan],los.ylos[ichan]]

  for i=0,nnpa-1 do oplot,[npaipos[i,0]],[npaipos[i,1]],psym=3,color=254./nnpa*i

  oplot, npafpos[*,0],npafpos[*,1],psym=3,thick=2

  ;; PLOT Histogram
  if keyword_set(ps) then device, filename='PLOTS/npa_histogram.eps' $
  else window,1

  ;; Define arrays for histograms
  ;; energy array:
  nen  = 51.
  emax = 100.
  emin = 0.
  energy_arr = emin+dindgen(nen)/(nen-1.)*(emax-emin)
  dE   = energy_arr[1]-energy_arr[0]

  ;;pitch array:
  npitch = 101
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
  histo =fltarr(ndet,nen)
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
              *fidasim.fbm.pitch_sign_convention
        dummy=min(abs(energy_arr-energy),eindex)
        dummy=min(abs(pitch_arr-pitch),pindex)    
        histo[idet,eindex]=histo[idet,eindex]+npawght[i]
        distri[idet,eindex,pindex]=distri[idet,eindex,pindex]+npawght[i]
        pEnergy[i] = eindex
     endfor

     contour,distri[idet,*,*],energy_arr,pitch_arr $
             ,c_colors=indgen(20)*12,nlevels=20,/fill,yran=[-1,1] $
             ,ytit='Pitch',xtit='Energy [keV]'
     plot, energy_arr, histo,ytit='neutrals/s',xtit='Energy [keV]'
  endfor
  if keyword_set(ps) then device, /close
end 

