pro transp_fbeam,inputs,coords,denf
;; Subroutine of FIDASIM to load a fast-ion distribution function
;; IN: inputs, coords
;; OUT: denf

  print, 'reading fast ion distribution function from transp output'
  cdftest=findfile(inputs.cdf_file)
  print, '======================='
  if cdftest[0] eq '' then begin
     print, " * File not found: ", inputs.cdf_file
     stop
  endif
  ;; ------------------------------
  ;; --load from TRANSP cdf file:--
  ;; ------------------------------
  cdfid=NCDF_Open(inputs.cdf_file,/nowrite)
  ;; Retrieve signals
  ;; --------------------------------------------------------
  ncdf_varget, cdfid,'TRANSP_RUNID', runid
  ncdf_varget, cdfid,'TIME' , cdf_time    
  ncdf_varget, cdfid,'R2D'  , r2d         ; rposition of cells
  ncdf_varget, cdfid,'Z2D'  , z2d         ;zposition of cells
  ncdf_varget, cdfid,'BMVOL', bmvol       ; plasma volume
  ncdf_varget, cdfid,'E_D_NBI', energy    ; central values
  ncdf_varget, cdfid,'A_D_NBI', pitch     ; central values
  ncdf_varget, cdfid,'F_D_NBI', FBM       ; fast-ion distribution function
  ncdf_varget, cdfid,'NTOT_D_NBI',ntot    ; total number of fast ions
  ncdf_varget, cdfid,'RSURF', rsurf       ; flux surface
  ncdf_varget, cdfid,'ZSURF', zsurf       ; flux surface
  NCDF_Close,cdfid
  ngrid=n_elements(r2d)
  ;;================================
  ;; get tranpped -passing boundary
  ;;================================
  rmin=fltarr(ngrid)
  for i=0,ngrid -1 do begin
     dummy=min((rsurf-r2d[i])^2+(zsurf-z2d[i])^2,index)
     index = array_indices(rsurf, index)
     rmin[i]=min(rsurf[*,index[1]])
  endfor
  pitch_boundary=sqrt(1.-rmin[*]/r2d[*])
  ;; ----------- Check the time
  print, ' CDF file time:',cdf_time 
  if abs(inputs.time-cdf_time) gt 0.02 then begin
     print, 'Time of CDF file and simulation disagree!'
     print, 'Time of CDF file and simulation disagree!'
     print, 'Time of CDF file and simulation disagree!'
     print, 'Time of CDF file and simulation disagree!'
     print, 'Time of CDF file and simulation disagree!'
                                ;print, '".c" to proceed anyway'
                                ;stop    
  endif     
  ;;-------------Convert eV-> keV
  energy=energy*1.0d-3  ;; fidasim needs energy in kev  
  fbm=fbm*1.0d3         ;; now, this needs to be corrected
  ;; as we now calculate with fast-ions/omega/keV/cm^3
  ;;------------Convert d_omega --> pitch
  ;; Fast-ion distribution is given as a function of cm^3, energy
  ;; and d_omega/4Pi. omega is the solild angle in 3D velocity space. In
  ;; order to transform this to a function depending on pitch instead
  ;; of d_omega/4PI, one has to multiply by 0.5!
  fbm*=0.5  
  ;; make sure that fbm is >=0:
  fbm>=0.
  ;;loading finished



  ;; TRANSP defines the pitch along the current direction. In
  ;; contrast, FIDASIM uses pitch along the B-field! Therefore,
  ;; reverse the pitch coordinate in fbm!
  npitch=n_elements(pitch)
  index=npitch-(indgen(npitch)+1)
  fbm[*,*,*]=fbm[*,index,*]
  ;;----------- select energy range -------
  index=where(energy ge inputs.emin and energy le inputs.emax,nenergy)
  energy=energy[index]     
  fbm=fbm[index,*,*]
  dE      = energy[2] - energy[1]
  emin=(float(energy[0])         - float(0.5*dE))>0.
  emax=float(energy[nenergy-1]) + float(0.5*dE)
  print, 'Energy min/max:', emin,emax
  ;; --------- select Pitch range --------
  index=where(pitch ge inputs.pmin and pitch le inputs.pmax,npitch)
  pitch=pitch[index] 
  fbm=fbm[*,index,*]   
  dP  = abs(pitch[2]  - pitch[1])
  pmin=(float(pitch[0])       - float(0.5*dP))>(-1)
  pmax=(float(pitch[npitch-1])+ float(0.5*dP))<1
  print, 'Pitch min/max:', pmin,pmax
  ;;-----------store distribution funciton in binary file
  file =inputs.root_dir+'RESULTS/'+inputs.fidasim_runid+'/transp_fbm.bin'
  openw, lun, file, /get_lun
  writeu,lun , strmid(inputs.cdf_file,strlen(inputs.cdf_file)-17,17)
  writeu,lun , double(cdf_time)
  ;; SPATIAL GRID
  writeu,lun , long(ngrid)
  writeu,lun , double(r2d[*])
  writeu,lun , double(z2d[*])
  writeu,lun , double(bmvol[*])
  ;; ENERGY GRID
  writeu,lun , long(nenergy)
  writeu,lun , double(emin)
  writeu,lun , double(emax)
  writeu,lun , double(energy[*])
  ;; PITCH GRID
  writeu,lun , long(npitch)
  writeu,lun , double(pmin)
  writeu,lun , double(pmax)
  writeu,lun , double(pitch[*])
  writeu,lun , double(FBM[*,*,*]) 
  close,lun
  free_lun, lun 
  print, 'TRANSP distribution in: '+file



  ;; ------map fdens on FIDASIM grid and sort out
  ;; ------points outside the separatrix
  ;;----------Determine fast-ion density averaged over pitch and energy
  fdens=total(reform(total(fbm,1)),1)*dE*dP
  denf=dblarr(coords.nx,coords.ny,coords.nz)*0.d0
  if ngrid le 220 then width=6. else width=4.
  for j=0,coords.ny-1 do begin
     rout=reform(coords.rrc_grid[*,j,0])
     zout=coords.zzc[*]
     TRIANGULATE, r2d, z2d, tr     
     fdens2=griddata(r2d,z2d,fdens,xout=rout,yout=zout $
                     ,/grid,/SHEPARDS,triangles=tr)
     for i=0L,coords.nx-1 do begin
        for k=0L,coords.nz-1 do begin
           a=sqrt((z2d-zout[k])^2+(r2d-rout[i])^2)
           ;; only write fdens2 if it is close to r2d,z2d grid
           if min(a) le width then denf[i,j,k]=fdens2[i,k] >0.        
        endfor
     endfor 
  endfor
end
