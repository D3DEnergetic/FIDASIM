pro start_fidasim  
  ;Definition of the simulation grid in machine coordinates
  nx=40
  ny=60
  nz=50
  xdim1=100  
  xdim2=220 
  ydim1=-190 
  ydim2=-10 
  zdim1=-75 
  zdim2=75 
  ;Basic grid points
  dx= (xdim2-xdim1) / double(nx) 
  dy= (ydim2-ydim1) / double(ny) 
  dz= (zdim2-zdim1) / double(nz)
  ;;cell borders
  xx=dx*dindgen(nx)+xdim1
  yy=dy*dindgen(ny)+ydim1
  zz=dz*dindgen(nz)+zdim1  


  user=+GETENV('USER')
  root_dir='/u/'+user+'/FIDASIM/'


  !Path = !path+ ':'+root_dir+'AUGD/'
  !Path = !path+ ':'+root_dir
  


  inputs={ doplot:[1], ps:[0] $ ;; plot, plot inputs into an eps file
           ,shot:28746L,time:4.42  $
           ,equil_exp:'AUGD' $
           ,equil_diag:'EQH' $
           ,rhostr:'rho_tor' $ ;; use input profiles als rho_pol or rho_tor
           ,root_dir:root_dir $
           ,fidasim_runid:'28746A01' $ ;; runid of F90FIDASIM
           ,isource:2 $   ;; Beam source (FIDASIM only simulates on NBI source)
           ,fida_diag:'CFR' $ ;; name of the FIDA diag
           ,npa:[0] $ ;; do a simulation for NPA

           ,no_spectra:[0] $ ;; if 1 then no spectra are calculated
           ,nofida:[0] $  ;; if 1 then no fast-ions are simulated
                       ;; Fast-ion distribution function from transp
           ,cdf_file:'TEST/28746A01_fi_1.cdf' $ ;; CDF file from transp with the distribution funciton
           ,emin:0. $  ;; minimum energy used from the distribution function
           ,emax:100. $;; maximum energy used from the distribution function
           ,load_neutrals:[0] $ ;; if 1 then the neutral density (prefida) is loaded from an existing neutrals.bin file
                       ;; dimensions for the simulation grid
           ,dx:dx,dy:dy,dz:dz,nx:nx,ny:ny,nz:nz,xx:xx,yy:yy,zz:zz $
           ,rotate:[0] $
           ,xdim1:xdim1 $  
           ,xdim2:xdim2 $  
           ,ydim1:ydim1 $
           ,ydim2:ydim2 $
           ,zdim1:zdim1 $
           ,zdim2:zdim2 $
                      ;; Number of Monte Carlo particles
           ,nr_fida:50000 $  ;; FIDA
           ,nr_ndmc:1000 $  ;; Beam emission
           ,nr_halo:50000 $  ;; Halo contribution
                      ;; Calculation of the weight function:
           ,calc_wght:[1] $ ;; if 1 then weight functions are calculated
           ,nr_wght:40 $  ;; Nr. of Pitches, energyies and gyro angles 
           ,emax_wght:100 $ ;; maximum energy
           ,ichan_wght:-1 $ ;; -1 for all channels, otherwise a given channel
           ,dwav_wght:1. $  ;; wavelength intervl
           ,wavel_start_wght:651. $ ;; minimum wavelength
           ,wavel_end_wght:663. }   ;; maximum wavelength
  

  ;; start the routine to collect the inputs
  dalpha_inputs,inputs
  ;; start the f90 routine
  result_dir = inputs.root_dir + 'RESULTS/' + inputs.fidasim_runid
  print, result_dir
;  spawn, 'module load intel'
;  spawn, root_dir+'fidasim '+result_dir


  ;; plot the spectra
;  plot_fidasim_spectra
  ;; plot the neutral densities
;  plot_fidasim_neutrals,/loga

end


