pro start_npasim  
  
  shot=28746
  time=4.421


  user=+GETENV('USER')
  root_dir='/u/'+user+'/FIDASIM2/'

  !Path = !path+ ':'+root_dir+'AUGD/'
  !Path = !path+ ':'+root_dir+'LIB/'
  !Path = !path+ ':'+root_dir
  

  npa_setup,det

  
  ;Definition of the simulation grid in machine coordinates
  nx=40
  ny=60
  nz=50
  xdim1=100.  
  xdim2=220. 
  ydim1=-150. 
  ydim2= 40. 
  zdim1=-75. 
  zdim2=75. 
  ;Basic grid points
  dx= (xdim2-xdim1) / double(nx) 
  dy= (ydim2-ydim1) / double(ny) 
  dz= (zdim2-zdim1) / double(nz)
  ;;cell borders
  xx=dx*dindgen(nx)+xdim1
  yy=dy*dindgen(ny)+ydim1
  zz=dz*dindgen(nz)+zdim1  
  ;; ------------------------
  ;; Define global variables:
  ;; ------------------------
  inputs={ doplot:[1], ps:[0] $   ;; plot, plot inputs into an eps file
           ,shot:shot,time:time  $
           ,equil_exp:'AUGD',equil_diag:'EQH',equil_ed:0 $
           ,rhostr:'rho_tor' $    ;; use input profiles als rho_pol or rho_tor
           ,root_dir:root_dir $
           ,fidasim_runid:string(shot,f='(i5)')+'npasim' $ ;;runid of F90FIDASIM
           ,isource:[2] $     ;; Beam source 
                      ;; Number of Monte Carlo particles
           ,nr_ndmc:5000 $  ;; Beam emission
           ,nr_halo:50000 $ ;; Halo contribution
            ;; dimensions for the simulation grid
           ,dx:dx,dy:dy,dz:dz,nx:nx,ny:ny,nz:nz,xx:xx,yy:yy,zz:zz $
           ,rotate:[0] $
           ,xdim1:xdim1 ,xdim2:xdim2 $  
           ,ydim1:ydim1 ,ydim2:ydim2 $
           ,zdim1:zdim1 ,zdim2:zdim2 $
           ;; Diagnostic settings
           ,diag:'NPA' $     ;; name of the FIDA diag
           ,det:det $
           ,npa:[1] $           ;; do a simulation for NPA
           ,calc_spec:[0] $     ;; if 1 then no spectra are calculated
           ,cdf_file:'TEST/28746A01_fi_1.cdf' $
           ,load_neutrals:[0] $ ;; load density from  existing file
           ,calc_birth: [0]   $ ;; calculate birth profile
           ;; FIDA/NPA simulation
           ,nr_fast:100000   $ ;; FIDA
           ,simfa:[1]        $ ;; simualte fast-ions
           ;; Fast-ion distribution function
           ,emin:0.   $ ;; minimum energy in FBM
           ,emax:100. $ ;; maximum energy in FBM
           ,pmin:-1.  $ ;; minimum pitch in FBM
           ,pmax:1.   $ ;; maximum pitch in FBM
           ;; weight functions:
           ,calc_wght:[0]    $ ;; if 1 then weight functions are calculated
           ,nr_wght:50       $ ;; Nr. of Pitches, energyies and gyro angles 
           ,ichan_wght:-1 $  ;; -1 for all channels, otherwise a given channel
           ,emax_wght:100. $ ;; maximum energy of weights
           ,dwav_wght:1. $   ;; wavelength interval
           ,wavel_start_wght:650. $ ;; min wavel of weight function
           ,wavel_end_wght:  665  }   ;; max wavel of weight function
  loadct,39 &  set_plot,'X' &  device, decomposed=0
  if inputs.ps[0] eq 1 then set_plot, 'ps'
  ;; -----------------------------------------------------
  ;; ------------- LOAD KINETIC PROFILES FOR FIDASIM -----
  ;; -----------------------------------------------------
  ;load_transp_profiles,inputs.shot,inputs.time,iprof,inputs.cdf_file
  ;save,filename='TEST/load_transp_profiles.idl',iprof
  restore,'TEST/load_transp_profiles.idl'
  ;; ----------------------------------------------------------
  ;; Start the IDL routine which collects and stores the inputs
  ;; ----------------------------------------------------------
  dalpha_inputs,inputs,iprof
  ;; start the f90 routine
  result_dir = inputs.root_dir + 'RESULTS/' + inputs.fidasim_runid
  spawn, root_dir+'fidasim '+result_dir

  ;; plot the spectra
  plot_npa,path=result_dir
end


