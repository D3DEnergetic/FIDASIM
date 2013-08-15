;; ------------------------------------------------------
;; ----------------- DALPHA -INPUTS ---------------------
;; ------------------------------------------------------
;; The main IDL routine of FIDASIM. It collects the inputs such as the
;; NBI geometry and power, kinetic profiles, equilibrium, diagnostic
;; geometry, fast-ion distribution function and stores the collected
;; data in a new FOLDER in the /RESULTS directory. THe file inputs.dat
;; and some binary files are written that are then read out by the
;; FORTRAN part of FIDASIM.
;; written by Benedikt Geiger 2013
pro los_track,coords,xyz_los_vec,xyspt,tcell,cell,ncell
  ri=xyspt
  vn=xyz_los_vec
  ;; zeros are not good!
  index=where(vn eq 0,nind)
  if nind gt 0 then vn[index]=0.0001
  p=[0,0,0]
  dummy = min( abs( coords.xxc - ri[0] ), index )
  p[0]=index
  dummy = min( abs( coords.yyc - ri[1] ), index ) 
  p[1]=index
  dummy = min( abs( coords.zzc - ri[2] ), index ) 
  p[2]=index
  tcellh=fltarr(1000)
  cellh=fltarr(3,1000)
  m=0
  cellh[*,m]=p    
;;loop along line of sight
  while(m lt 1000) do begin
     l=p
     index=where(vn gt 0.d0)
     if index[0] ne -1 then begin
        l(index)=p(index)+1
     endif 
     if l[0] gt coords.nx-1 or l[1] gt coords.ny-1 or $
        l[2] gt coords.nz-1 then goto, out
     ;;time needed to go into next cell
     dt_arr=fltarr(3)
     dt_arr[0] = ( coords.xx(l[0]) - ri[0] ) /vn[0]
     dt_arr[1] = ( coords.yy(l[1]) - ri[1] ) /vn[1]
     dt_arr[2] = ( coords.zz(l[2]) - ri[2] ) /vn[2]     
     dt=min(dt_arr,index)
     ri[*] = ri[*] + vn[*]*dt 
     if vn[index] gt 0.d0 then begin
        p[index]=p[index]+1
     endif else begin
        p[index]=p[index]-1
     endelse
     if p[0] lt 0 or p[1] lt 0 or p[2] lt 0 then goto, out
     tcellh[m]=dt
     m=m+1
     cellh[*,m]=p[*] 
  end  
  out:
;; Store results into compressed arrays!
  if m gt 1 then begin
     tcell  = tcellh[0:m-2]
     cell   = cellh[*,0:m-2]
  endif else begin
     tcell   = -1
     cell   = -1
  endelse
  ncell=m-1
end



pro dalpha_inputs,inputs,profiles
;Original version from W. W. Heidbrink 2007
;rewritten by Benedikt Geiger 2012 (IPP-Garching, ASDEX Upgrade)
;routine to collect and store inputs for the fortran routine (fidasim2.f90)
;IN: inputs, profiles
  compile_opt defint32
  ;;----------------------------------------------------------
  ;;GEOMETRY/DATA, FIDAsim grid, NUMERICS, WAVELENGTH grid
  ;;----------------------------------------------------------
  BTIPsign=-1.d0              ; Bt and Ip are in the opposite direction   
  ab=2.01410178d0             ; atomic mass of beam [u]
  ai=2.01410178d0             ; atomic mass of hydrogenic plasma ions [u]
  impurity_charge=6           ; 5: BORON, 6: carbon, 7: Nitrogen               
  ;;------------------------------------------------
  ;; Wavelength grid  
  lambdamin=6470.d0           ; minimum wavelength of wavelength grid[A] 
  lambdamax=6670.d0           ; maximum wavelength of wavelength grid[A] 
  nlambda=2000L               ; number of wavelengths
  dlambda= (lambdamax-lambdamin)/float(nlambda)                            
  ;;---------------------------------------------------------------------
  ;; FIDASIM grid   
  nx=long(inputs.nx) & ny=long(inputs.ny) & nz=long(inputs.nz)
  ng=nx*ny*nz     ;; nr of cells
  dx=inputs.dx & dy=inputs.dy & dz=inputs.dz
  dr=[dx,dy,dz]        ;; size of cells
  drmin=min(dr)        ;; minimal size
  dv=dr[0]*dr[1]*dr[2] ;;volume
  ;; Basic grid points
  ;; cell borders
  xx=inputs.xx & yy=inputs.yy & zz=inputs.zz
  ;; cell centers
  xxc=xx+0.5d0*dx & yyc=yy+0.5d0*dy & zzc=zz+0.5d0*dz
  ;; Major radius and toroidal angle
  rrc_grid=dblarr(nx,ny,nz)
  phi_grid=dblarr(nx,ny,nz)
  zzc_grid=dblarr(nx,ny,nz)
  for i=0L,nx-1 do begin
     for j=0L,ny-1 do begin
        for k=0L,nz-1 do begin
           rrc_grid[i,j,k]=sqrt(xxc[i]^2+yyc[j]^2)
  
           phi_grid[i,j,k]=atan(yyc[j]/xxc[i])
           if xxc[i] gt 0 and yyc[j] lt 0 then begin
              phi_grid[i,j,k]+=2.*!pi
           endif else begin
              if xxc[i] lt 0 and yyc[j] gt 0 then begin
              phi_grid[i,j,k]+=!pi
              endif else begin
                 if xxc[i] lt 0 and yyc[j] lt 0 then phi_grid[i,j,k]+=!pi
              endelse
           endelse
           zzc_grid[i,j,k]=zzc[k]
        endfor
     endfor
  endfor
  ;;----------------------------------------------------------
  ;; NBI GEOMETRY and DATA
  ;nbi_data,inputs,einj,pinj $ ;power,energy,particle mix
  ;              ,ffull,fhalf,fthird,doplot=0 $
  ;              ,ps=inputs.ps  
  ;save,filename='TEST/nbi_data.idl',einj,pinj,ffull,fhalf,fthird
  restore,'TEST/nbi_data.idl'
  ;; check if selected source is on!
  if pinj[inputs.isource] le 0. then begin
     print, 'the selected source nr',inputs.isource,' is not on!'
     stop
  endif

  ;;fetch TRANSP NBI geometry
  nbi_geometry_transp,nbgeom,rotate=inputs.rotate
  ;;------------------------------------------------------------------
  ;; - STORE GEOMETRY AND GENERAL INPUTS IN STURCTURE AND "inputs.dat"
  ;;------------------------------------------------------------------
  coords={dx:dx,dy:dy,dz:dz,drmin:drmin,dv:dv,ng:ng,nx:nx,ny:ny,nz:nz $
          ,xx:xx,yy:yy,zz:zz,xxc:xxc,yyc:yyc, zzc:zzc $
          ,rrc_grid:rrc_grid, phi_grid:phi_grid, zzc_grid:zzc_grid $
          ,isource:inputs.isource}

  runid=inputs.fidasim_runid
  if file_test(inputs.root_dir+'RESULTS/'+runid,/directory) eq 0 then begin
     spawn,'mkdir '+inputs.root_dir+'RESULTS/'+runid
  endif

  file = inputs.root_dir+'RESULTS/'+runid+'/inputs.dat'
  openw, 55, file
  printf,55,'# FIDASIM input file created: ', systime(),' Version 2.1'
  printf,55, inputs.root_dir
  printf,55, inputs.shot         ,f='(i6,"         # shotnumber")'  
  printf,55, inputs.time,f='(1f8.5,"       # time")'
  printf,55, runid
  printf,55,' ',inputs.diag, '           # diagnostic'
  printf,55,inputs.calc_birth,f='(i2,"             # calculate birth profile")'
  printf,55,inputs.calc_spec,f='(i2,"             # calculate spectra")'
  printf,55,inputs.ps,f='(i2,"             # plot ps output")'
  printf,55,inputs.npa          ,f='(i2,"             # NPA simulation")'
  printf,55,inputs.load_neutrals,f='(i2,"             # load NBI+HALO density")'
  printf,55, inputs.guidingcenter,f='(i2,"             # 0 for full-orbit F")'
  printf,55, inputs.f90brems,f='(i2,"             # 0 reads IDL v.b.")'
  printf,55,inputs.calc_wght,f='(i2,"             # calculate wght function")'
  printf,55,'# weight function settings:'
  printf,55,inputs.nr_wght,f='(i9,"      # number velocities")'
  printf,55,inputs.ichan_wght,f='(i3,"      # channel for weight function")'
  printf,55,inputs.emax_wght,f='(1f12.2,"       # emax for weights")'
  printf,55,inputs.dwav_wght,f='(1f12.5,"       # dwav")'
  printf,55,inputs.wavel_start_wght,f='(1f12.5,"       # wavel_start")'
  printf,55,inputs.wavel_end_wght,f='(1f12.5,"       # wavel_end")'
  printf,55,'# Monte Carlo settings:'
  printf,55,inputs.nr_fast,f='(i9,"      # number of FIDA mc particles")'  
  printf,55,inputs.nr_ndmc,f='(i9,"      # number of NBI mc particles")' 
  printf,55,inputs.nr_halo,f='(i9,"      # number of HALO mc particles")'
  printf,55,impurity_charge,f='(i2,"             # Impurity charge")'
  printf,55,'# Location of transp cdf file:'
  printf,55,inputs.cdf_file
  printf,55,'# discharge parameters:'
  printf,55,btipsign,f='(i3,"            # B*Ip sign")'
  printf,55,ai,f='(1f7.4,"        # plasma mass")'
  printf,55,ab,f='(1f7.4,"        # NBI mass")'
  printf,55,'# wavelength grid:'
  printf,55,nlambda,f='(1i5,"          # nlambda")'
  printf,55,lambdamin,f='(1f9.3,"      # lambda min")'
  printf,55,lambdamax,f='(1f9.3,"      # lambda max")'
  printf,55,'# simulation grid: '
  printf,55,nx,f='(1i3,"            # nx")'
  printf,55,ny,f='(1i3,"            # ny")'
  printf,55,nz,f='(1i3,"            # nz")'  
  for i=0L,nx-1 do begin   ;; cell borders          
     printf,55,xx[i],f='(1f9.4,"      # xx[i]")'
  endfor
  for i=0L,ny-1 do begin
     printf,55,yy[i],f='(1f9.4,"      # yy[i]")'
  endfor
  for i=0L,nz-1 do begin
     printf,55,zz[i],f='(1f9.4,"      # zz[i]")'
  endfor
  printf,55,'# Neutral beam injection:'
  printf,55,nbgeom.BMWIDRA,f='(1f9.4,"      # NBI half width horizontal")'
  printf,55,nbgeom.BMWIDZA,f='(1f9.4,"      # NBI half width vertical")'
  ii=inputs.isource[0]
  printf,55, ii,f='(1i2,"             # Nr of NBI")'
  printf,55,nbgeom.divy[0,ii],f='(1f10.7,"     #divergence y of full comp")'
  printf,55,nbgeom.divy[1,ii],f='(1f10.7,"     #divergence y of half comp")'
  printf,55,nbgeom.divy[2,ii],f='(1f10.7,"     #divergence y of third comp")'
  printf,55,nbgeom.divz[0,ii],f='(1f10.7,"     #divergence z of full comp")'
  printf,55,nbgeom.divz[1,ii],f='(1f10.7,"     #divergence z of half comp")'
  printf,55,nbgeom.divz[2,ii],f='(1f10.7,"     #divergence z of third comp")'
  printf,55,nbgeom.focy[ii],f='(1f9.4,"      # focal length in y")' 
  printf,55,nbgeom.focz[ii],f='(1f9.4,"      # focal length in z")' 
  printf,55,einj[ii],f='(1f9.4,"      # injected energy [keV]")' 
  printf,55,pinj[ii],f='(1f9.4,"      # injected power [MW]")'  
  printf,55,'# Species-mix (Particles):'
  printf,55 ,ffull[ii],f='(1f9.6,"      # full energy")' 
  printf,55 ,fhalf[ii],f='(1f9.6,"      # half energy")'  
  printf,55 ,fthird[ii],f='(1f9.6,"      # third energy")' 
  printf,55, '#position of NBI source in xyz coords:'
  printf,55,nbgeom.xyz_src[ii,0],f='(1f9.4,"      # x [cm]")' 
  printf,55,nbgeom.xyz_src[ii,1],f='(1f9.4,"      # y [cm]")' 
  printf,55,nbgeom.xyz_src[ii,2],f='(1f9.4,"      # z [cm]")' 
  printf,55,'# 3 rotation matrizes 3x3'
  for j=0,2 do begin
     for k=0,2 do begin
        printf,55 ,nbgeom.Arot[ii,j,k] ;; rotation in the top-down view plane
        printf,55 ,nbgeom.Brot[ii,j,k] ;; vertical rotation
        printf,55 ,nbgeom.Crot[ii,j,k] ;; vertical rotation
     endfor 
  endfor
  printf,55, '#rho-string:'
  printf,55,profiles.rho_str
  printf,55, '#rotate: (for different coordinate system)'
  printf,55,inputs.rotate,f='(1f10.5,"     #rotate")'
  close,55


  ;;-------------------------------------------------
  ;; MAP kinetic profiles on FIDASIM grid
  ;;------------------------------------------------- 
  ;; print, 'time:' ,inputs.time    
  ;; Magnetic field (bx,by,bz) and RHO
 ; load_bfield,inputs,coords,bx,by,bz,rhopf,rhotf
 ; if profiles.rho_str eq 'rho_tor' then rho_grid=rhotf 
 ; if profiles.rho_str eq 'rho_pol' then rho_grid=rhopf
 ; save,filename='TEST/load_bfield.idl',rho_grid,bx,by,bz
  restore,'TEST/load_bfield.idl'
  ;;Electric field
  ex=bx*0.0d0  
  ey=by*0.0d0
  ez=bz*0.0d0
  ;;Electron density
  dene      = 1.d-6 * interpol(profiles.dene,profiles.rho,rho_grid)>0. ;[1/cm^3]
  ;;Zeff
  zeff      = interpol(profiles.zeff,profiles.rho,rho_grid) >1.
  ;;Impurity density
  deni = (zeff-1.)/(impurity_charge*(impurity_charge-1))*dene
  ;;Proton density
  denp = dene-impurity_charge*deni
  print,total(deni)/total(denp)*100. ,'percent of impurities'
  ;;Fast-ion density
  if inputs.simfa eq 1 and inputs.nr_fast gt 1 then begin 
     transp_fbeam,inputs,coords,denf
  endif else begin
     denf=dene*0.d0
  endelse
  ;;Electron temperature
  te       = 1.d-3 * interpol(profiles.te,profiles.rho,rho_grid) ;keV
  if min(te) lt 0 then te[where(te lt 0.)]=0.d0
  ;;Ion temperature   
  ti       = 1.d-3 * interpol(profiles.ti,profiles.rho,rho_grid)>0. ;keV
  if max(ti) gt 20. or max(te) gt 20. then begin
     print, 'Look at the tables, they might only consider'
     print, 'temperatures unitl 10keV!'
     stop
  endif
  ;;Plasma rotation
  vtor     = 1.d2* interpol(profiles.vtor,profiles.rho,rho_grid) ; [cm/s]  
  vrotx = - cos(!pi*0.5d0-coords.phi_grid)*vtor 
  vroty =   sin(!pi*0.5d0-coords.phi_grid)*vtor
  vrotz =   vrotx* 0.d0 


  ;; test if there are NANs or Infinites in the input profiels
  index=where(finite([ti,te,dene,denp,zeff,denp,deni]) eq 0,nind)
  if nind gt 0 then stop


  ;;-------SAVE-------
  plasma={rho_grid:rho_grid $
          ,bx:bx, by:by, bz:bz $
          ,ab:ab,ai:ai,einj:einj,pinj:pinj $
          ,te:te,ti:ti,vtor:vtor,dene:dene,denp:denp $
          ,deni:deni,denf:denf,zeff:zeff}  
  

  file =inputs.root_dir+'RESULTS/'+runid+'/plasma.bin'
  openw, lun, file, /get_lun
  writeu,lun , long(coords.nx)
  writeu,lun , long(coords.ny)
  writeu,lun , long(coords.nz)
  writeu,lun , double(te)
  writeu,lun , double(ti)
  writeu,lun , double(dene)
  writeu,lun , double(denp)
  writeu,lun , double(deni)
  writeu,lun , double(denf)
  writeu,lun , double(vrotx)
  writeu,lun , double(vroty)
  writeu,lun , double(vrotz)
  writeu,lun , double(zeff)
  writeu,lun , double(bx)
  writeu,lun , double(by)
  writeu,lun , double(bz)
  writeu,lun , double(ex)
  writeu,lun , double(ey)
  writeu,lun , double(ez)
  writeu,lun , double(rho_grid)
  close,lun
  free_lun, lun
  ;;print, 'plasma parameters stored in BINARY: '+file  

  ;;-----------------------------------------------------------
  ;;FIDA diagnositc ----- Detector vectors, weights and spectra
  ;;-----------------------------------------------------------
  if inputs.calc_spec or inputs.calc_wght or inputs.npa then begin
     det=inputs.det
     ;; the sigma_pi ratio is the ratio between the pi and sigma lines
     ;; caused by the stark broadening. If there are e.g. mirrors used
     ;; in the optics, then this ratio can change! 
     ;; the sigma_pi_ratio is multiplied with the sigma component
     ;; before normalizing the intensity of the stark components (see
    ;; function spectrum in fidasim.f90)
    ;;CALCULATE WEIGHTS
     weight  = replicate(0.d0,coords.nx,coords.ny,coords.nz,det.nchan)
     print, 'nchan:', det.nchan
     for chan=0, det.nchan-1 do  begin
        xyzhead = [det.xhead[chan],det.yhead[chan],det.zhead[chan]]
        xyzlos  = [det.xlos[chan], det.ylos[chan], det.zlos[chan]]
        vi    = xyzlos-xyzhead
        dummy = max(abs(vi),ic)
        nstep = fix(700./dr[ic])
        vi    = vi/sqrt(vi[0]^2+vi[1]^2+vi[2]^2) ;; unit vector
        xyz_pos = xyzhead
      ; find first grid cell
        for i=0,nstep do begin
           xyz_pos[0] = xyz_pos[0] + dr[ic] * vi[0]/abs(vi[ic])
           xyz_pos[1] = xyz_pos[1] + dr[ic] * vi[1]/abs(vi[ic])
           xyz_pos[2] = xyz_pos[2] + dr[ic] * vi[2]/abs(vi[ic])
           if xyz_pos[0] gt xx[0] and xyz_pos[0] lt xx[nx-1]+dx and $ 
              xyz_pos[1] gt yy[0] and xyz_pos[1] lt yy[ny-1]+dy and $
              xyz_pos[2] gt zz[0] and xyz_pos[2] lt zz[nz-1]+dz then begin
              goto, out
           endif  
        endfor
        out:
      ; determine cells along the LOS
        if i lt nstep then begin
           los_track,coords,vi,xyz_pos,tcell,cell,ncell
           if ncell gt 1 then begin
              for jj=0,ncell-1 do begin
                 if finite(tcell[jj]) eq 0 then stop
                 ;;  tcell is the length of the track (cm) as v is 1cm/s
                 weight[cell[0,jj],cell[1,jj],cell[2,jj],chan]=tcell[jj]
              endfor
           endif else begin
              print, 'LOS only crosses one cell!'
           endelse
        endif else begin
           print, 'LOS does not cross the simulation grid!'
           ;;print,'chan: ', chan
        endelse
     endfor
     index=where(finite(weight) eq 0,nind)
     if nind gt 0 then begin
        print,'weight set to 0. as it was NAN or Infinite!'
        weight[index]=0.
     endif
     detector={ det:det  , weight:weight }
  endif else begin 
     det=set_det()
     detector={det:det,weight:replicate(0.d0,coords.nx,coords.ny,coords.nz,1) }
  endelse
  ;;-------SAVE-------
  file =inputs.root_dir+'RESULTS/'+runid+'/los.bin'
  openw, lun, file, /get_lun
  writeu,lun , long(det.nchan)
  for chan=0,det.nchan-1 do begin
     writeu,lun, double(det.xhead[chan])
     writeu,lun, double(det.yhead[chan])
     writeu,lun, double(det.zhead[chan])
     writeu,lun, double(det.xlos[chan])
     writeu,lun, double(det.ylos[chan])
     writeu,lun, double(det.zlos[chan])
     writeu,lun, double(det.npa.headsize[chan])
     writeu,lun, double(det.npa.opening_angle[chan])
     writeu,lun, double(det.sigma_pi[chan])
  endfor
  writeu,lun , double(detector.weight[*,*,*,*])
  close,lun
  free_lun, lun
  ;print, 'LOS parameters stored in BINARY: '+file 
  

  ;;--------------
  ;;Plot inputs
  ;;--------------
  if inputs.doplot then plot_dalpha_inputs,inputs,nbgeom,coords,plasma,detector,dens=dens
  ; print, 'inputs stored!',inputs.shot,inputs.time
  ; print,systime()
  if inputs.ps then device,/close
end
 



