pro transp_fbeam,inputs,coords,denf
  !P.charsize=1.
  !P.background=255 & !P.color=0
  doplot=1
  print, 'reading fast ion distribution function from transp output'
  cdftest=findfile(inputs.cdf_file)
  print, '======================='
  if cdftest[0] ne '' then  print, " * Found file: ", inputs.cdf_file
  cdfid=NCDF_Open(inputs.cdf_file,/nowrite)
  ;; Retrieve signals
  ;; --------------------------------------------------------
  ncdf_varget, cdfid,'TRANSP_RUNID', runid
  ncdf_varget, cdfid,'TIME' , cdf_time    
  ncdf_varget, cdfid,'R2D'  , r2d     ; rposition of cells
  ncdf_varget, cdfid,'Z2D'  , z2d     ;zposition of cells
  ncdf_varget, cdfid,'BMVOL', bmvol   ; plasma volume
  ncdf_varget, cdfid,'E_D_NBI', energy ; central values
  ncdf_varget, cdfid,'A_D_NBI', pitch  ; central values
  ncdf_varget, cdfid,'F_D_NBI', FBM    ; fast-ion distribution function
  ncdf_varget, cdfid,'NTOT_D_NBI',ntot ; total number of fast ions
  ncdf_varget, cdfid,'RSURF', rsurf    ; flux surface
  ncdf_varget, cdfid,'ZSURF', zsurf    ; flux surface
  NCDF_Close,cdfid
  ;; ----------- Check the time
  print, ' CDF file time:',cdf_time 
  if abs(inputs.time-cdf_time) gt 0.02 then begin
     print, 'Time of CDF file and simulation disagree!'
     ;print, '".c" to proceed anyway'
     ;stop    
  endif     
  ;;-------------Convert eV-> keV
  energy=energy*1.0d-3          ;; fidasim needs energy in kev  
  fbm=fbm*1.0d3                 ;; now, this needs to be corrected
  ;; as we now calculate with fast-ions/omega/keV/cm^3
                                ;;------------Convert d_omega --> pitch
  ;; Fast-ion distribution is given as a function of cm^3, energy
  ;; and d_omega/4Pi. omega is the solild angle in 3D velocity space. In
  ;; order to transform this to a function depending on pitch instead
  ;; of d_omega/4PI, one has to multiply by 0.5!
  fbm=fbm*0.5
  ;;------------Cut off energy
  print,'EMIN: ',inputs.emin
  print,'EMAX: ',inputs.emax
  w=where(energy ge inputs.emin and energy le inputs.emax)
  energy=energy[w]     
  fbm=fbm[w,*,*]
  index=where(fbm lt 0.,nind) 
  if nind gt 0. then fbm[index]=0.d0





  ;;-----------store distribution funciton in binary file
  file ='RESULTS/'+inputs.fidasim_runid+'/transp_fbm.bin'
  sz=size(FBM)
  nenergy=sz(1)
  npitch=sz(2)
  ngrid=sz(3)
  openw, lun, file, /get_lun
  writeu,lun , long(ngrid)
  for i=0,ngrid-1 do writeu ,lun , float(r2d[i])
  for i=0,ngrid-1 do writeu ,lun , float(z2d[i])
  writeu,lun , long(nenergy)
  writeu,lun , double(inputs.emin)
  writeu,lun , double(inputs.emax)
  for i=0,nenergy-1 do writeu ,lun , double(energy[i])
  writeu,lun , long(npitch)
  writeu,lun , double(pitch[0])
  writeu,lun , double(pitch[npitch-1])
  for i=0,npitch-1 do writeu ,lun , double(pitch[i])
  for i=0,nenergy-1 do begin
     for j=0,npitch-1 do begin
        for k=0,ngrid-1 do begin    
           writeu ,lun , float(FBM[i,j,k]/max(FBM[*,*,k])) ;; normalized
        endfor
     endfor
  endfor
  close,lun
  free_lun, lun
  print, 'TRANSP distribution in: '+file




  ;;----------Determine fast-ion density averaged over pitch and energy
  dE      = energy[2] - energy[1]
  dpitch  = pitch[2]  - pitch[1]
  fdens=total(reform(total(fbm,1)),1)*dE*dpitch




  ;; ------map fdens on FIDASIM grid and sort out
  ;; ------points outside the separatrix
  ;;------Determine FIDAsim grid to map distribution function
  r_grid=fltarr(coords.nx,coords.ny,coords.nz)
  for i=0L,coords.nx-1 do begin
     for j=0L,coords.ny-1 do begin
        for k=0L,coords.nz-1 do begin
           jj=i+coords.nx*j+coords.nx*coords.ny*k
          r_grid[i,j,k]=coords.r_grid[jj]
       endfor      
     endfor
  endfor 
  f3d=dblarr(coords.nx,coords.ny,coords.nz)*0.d0
  if ngrid le 220 then width=6. else width=4.
  for j=0,coords.ny-1 do begin
     rout=reform(r_grid[*,j,0])
     zout=coords.zzc[*]
     TRIANGULATE, r2d, z2d, tr     
     fdens2=griddata(r2d,z2d,fdens,xout=rout,yout=zout $
                     ,/grid,/SHEPARDS,triangles=tr)
     for i=0L,coords.nx-1 do begin
        for k=0L,coords.nz-1 do begin
           r=rout[i]
           z=zout[k]
           a=sqrt((z2d-z)^2+(r2d-r)^2)
           amin=min(a)
           ;; only write fdens2 if it is close to r2d,z2d grid
           if amin le width then f3d[i,j,k]=fdens2[i,k]         
        endfor
     endfor 
  endfor
  denf=dblarr(coords.ng)*0.d0   
  for i=0L,coords.nx-1 do begin
     for j=0L,coords.ny-1 do begin
        for k=0L,coords.nz-1 do begin
           jj=i+coords.nx*j+coords.nx*coords.ny*k
           denf[jj]=f3d[i,j,k] >0.
       endfor      
     endfor
  endfor 


  
end
