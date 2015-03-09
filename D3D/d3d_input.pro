PRO d3d_input,inputs                                   

  ;;-----------------------------------------------------
  ;;				PREFIDA INPUT FILE
  ;;-----------------------------------------------------
  ;; Comment
  comment = 'This is a template for D3D runs'

  ;; LONG Shot Number
  shot=157545L

  ;; FLOAT Time (seconds) has to match GAProfiles time exactly.
  time=3.750

  ;; STRING runid of FIDASIM (10 characters or less)
  runid='def'

  ;; STRING D3D,NSTX,AUGD,MAST
  device='D3D'

  ;; STRING Name of the FIDA diag
  ;; Valid names are
  ;; MAIN_ION30
  ;; MAIN_ION210
  ;; MAIN_ION330
  ;; TANGENTIAL
  diag='MAIN_ION30'

  ;; STRING Location where results will be stored 
  ;; result_dir/runid directory will be created
  ;; Needs trailing '/'
  result_dir='/u/grierson/FIDASIM/RESULTS/D3D/'+$
             STRTRIM(shot,2)+'/'+$
             STRING(time*1.e3,FOR='(I05)')+'/'+$
             diag+'/'

  ;; STRING Location of profile save files.
  ;; Needs trailing '/'
  profile_dir='/u/grierson/gaprofiles/f90fidasim/'+$
              STRTRIM(shot,2)+'/'+$
              STRING(time*1.e3,FOR='(I05)')+'/'+$
              diag+'/'+$
              runid+'/'

  ;;----------------------------------------------------
  ;; Fast-ion distribution function from transp
  ;;----------------------------------------------------
  ;; STRING CDF file from TRANSP+NUMBEM with distribution function
  cdf_file='/u/grierson/transp/fbm/blank_fi_1.cdf'

  ;; FLOAT minimum energy used from the distribution function
  emin=10.

  ;; FLOAT maximum energy used from the distribution function
  emax=100.

  ;; FLOAT Minimum pitch used from the distribution function (depricated)
  pmin=-1.

  ;; FLOAT Maximum pitch used from the distribution function (depricated)
  pmax=1.

  ;;-----------------------------------------------------
  ;; Beam/diagnostic/equilibrium Selection
  ;;-----------------------------------------------------
  ;; INTEGER Beam source index (FIDASIM only simulates one NBI source)
  ;; 0 = 30LT
  ;; 1 = 30RT
  ;; 4 = 210LT
  ;; 5 = 210RT
  ;; 6 = 330LT
  ;; 7 = 330RT
  isource=6

  ;; FLOAT Accelerator voltage [keV].
  ;; If 0, get data from MDS+
  einj=75.26

  ;; FLOAT Injected Power [MW].
  ;; If 0, get data from MDS+
  pinj=1.915
  
  ;; GEQDSK file else ''
  ;; Depricated
  ;; Will look for GEQDSK in profiles directory first.
;  gfile=''

  ;; STRING MDS+ equilibrium runid
  equil='EFIT02'

  ;;-----------------------------------------------------
  ;; Discharge Parameters
  ;;-----------------------------------------------------
  ;; DOUBLE Sign of Bt*Ip
  ;; Use -1 when Bt and Ip are in the opposite direction     
  btipsign=1.d0 

  ;; DOUBLE Atomic mass of beam [u]
  ab=2.01410178d0 

  ;; DOUBLE Atomic mass of hydrogenic plasma ions [u]
  ai=2.01410178d0
  
  ;; INTEGER Impurity charge
  ;; 5: BORON, 6: carbon, 7: Nitrogen
  impurity_charge=6 

  ;;-----------------------------------------------------
  ;; Wavelength Grid
  ;;-----------------------------------------------------
  IF STRCMP(diag,'MAIN_ION30') OR $
    STRCMP(diag,'MAIN_ION210') OR $
    STRCMP(diag,'MAIN_ION330') THEN BEGIN   

      ;; Main-ion spectrometers tuned to D-alpha

      ;; DOUBLE Minimum wavelength of wavelength grid [nm]
      lambdamin=649.d0
      
      ;; DOUBLE Maximum wavelength of wavelength grid[nm] 
      lambdamax=663.d0
      
      ;; LONG Number of wavelengths
      nlambda=768L
  ENDIF

  IF STRCMP(diag,'ET330') THEN BEGIN   

      ;; ET systems tuned edge FIDA

      ;; DOUBLE Minimum wavelength of wavelength grid [nm]
      lambdamin=657.d0
      ;; DOUBLE Maximum wavelength of wavelength grid[nm] 
      lambdamax=663.d0
      ;; LONG Number of wavelengths
      nlambda=326L
  ENDIF

  ;; Wavelength separation
  dlambda= (lambdamax-lambdamin)/double(nlambda)	

  ;;---------------------------------------------------
  ;; Define FIDASIM grid in machine coordinates(x,y,z)
  ;;---------------------------------------------------
  ;; From FIDAcode_rel_3_2 alpha_NB
  ;; [-2.51013, -2.35899, 1.67866, 1.82980, 1.31179, 1.46293, -1.46293, -1.31179]
  IF isource EQ 0 OR isource EQ 1 THEN BEGIN
      ;; FLOAT 30LT/RT crossover (u,v,z) [cm]
      origin=[129.26, 239.35,0.0]

      ;; Rotation of beam about z-axis.
      ;; Positive angle is CCW from east wall (90 deg port)
      ;; 30LT
      IF isource EQ 0 THEN alpha = -2.5206D
      ;; This second one is the old IDL FIDAsim angle.
;      IF isource EQ 0 THEN alpha = -2.51013D
      ;; 30RT
      IF isource EQ 1 THEN alpha = -2.3590D
      
      ;; Beam tilt angle
      beta=0
  ENDIF

  IF isource EQ 6 OR isource EQ 7 THEN BEGIN
      ;; FLOAT 330LT/RT crossover (u,v,z) [cm]
      origin=[-139.49, 233.77, 0.0]
      
      ;; Rotation of beam about z-axis.
      ;; Positive angle is CCW from east wall (90 deg port)
      ;; 330LT
      IF isource EQ 6 THEN alpha = -1.46293D
      
      ;; 330RT alpha
      IF isource EQ 7 THEN alpha=-1.31179D

      beta=0
  ENDIF

  ;; Spatial domain and resolution
  ;; For 30 and 210 degree beams core use xmin=40, xmax=140
  ;; For 330 degree beam edge use xmin=20, xmax=80
  IF isource LT 6 THEN BEGIN
      ;; Number of cells in x direction      
      nx=48
      ;; Number of cells in y direction
      ny=48
      ;; Number of cells in z direction
      nz=50
      ;; Minimum x value
      xmin=40.
      ;; Maximum x value
      xmax=140.
      ;; Minimum y value
      ymin=-50.
      ;; Maximum y value
      ymax=50.
      ;; Minimum z value
      zmin=-70.
      ;; Maximum z value
      zmax=70.
  ENDIF
  IF isource GE 6 THEN BEGIN
      ;; Number of cells in x direction      
      nx=48*2
      ;; Number of cells in y direction
      ny=48
      ;; Number of cells in z direction
      nz=50
      ;; Minimum x value
      xmin=20.
      ;; Maximum x value
      xmax=80.
      ;; Minimum y value
      ymin=-50.	
      ;; Maximum y value
      ymax=50.
      ;; Minimum z value
      zmin=-70.
      ;; Maximum z value
      zmax=70.
  ENDIF


  ;;--------------------------------------------------
  ;; Define number of Monte Carlo particles
  ;;--------------------------------------------------
  ;; FIDA
  nr_fast=5000000
  ;; Beam particles
  nr_nbi=100000
  ;; Halo particles
  nr_halo=1000000

  
  ;;--------------------------------------------------
  ;; Calculation of the weight function
  ;;--------------------------------------------------
  ;; Number of energies
  ne_wght=50
  ;; Numbre of pitches
  np_wght=50
  ;; Number of gyro-angles
  nphi_wght=50
  ;; Maximum energy (keV)
  emax_wght=100.
  ;; -1 for all channels, otherwise a given channel index
  ichan_wght=-1
  ;; Wavelength interval
  dwav_wght=1.
  ;; Minimum wavelength
  wavel_start_wght=651.
  ;; Maximum wavelength
  wavel_end_wght=663.
  

  ;;-------------------------------------------------
  ;; Simulation switches
  ;; These are binary (0 or 1) to perform calcualtion or not.
  ;;-------------------------------------------------
  ;; NPA is calcuated?
  calc_npa=[0]

  ;; Spectra are calculated ?
  calc_spec=[1]

  ;; The birth profile is calculated?
  calc_birth=[0]

  ;; If 0 use the IDL bremstrahlung calculation rather than f90
  calc_brems=[0]

  ;; If 1 then weight functions are calculated
  calc_fida_wght=[0]

  ;; If 1 then weight functions are calculated
  calc_npa_wght=[0]

  ;; If 1 then the neutral density is loaded from an existing
  ;; neutrals.bin file located in runid directory
  load_neutrals=[0]

  ;; Load fbm
  load_fbm=[1]
  
  ;; If interactive then the % completion is shown
  interactive=[0]

  ;;------------------------------------------------
  ;; DO NOT MODIFY THIS PART
  ;;------------------------------------------------
  install_dir=+GETENV('FIDASIM_DIR')
  
  inputs={comment:comment,shot:shot,time:time,runid:runid,device:strupcase(device),install_dir:install_dir,result_dir:result_dir,$
          cdf_file:cdf_file,profile_dir:profile_dir,emin:emin,emax:emax,pmin:pmin,pmax:pmax,isource:isource,diag:diag,$
          einj:einj,pinj:pinj,equil:equil,btipsign:btipsign,ab:ab,ai:ai,impurity_charge:impurity_charge,$
          lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,$
          nx:nx,ny:ny,nz:nz,xmin:xmin,xmax:xmax,ymin:ymin,ymax:ymax,zmin:zmin,zmax:zmax,$
          origin:origin,alpha:alpha,beta:beta,$
          nr_fast:nr_fast,nr_nbi:nr_nbi,nr_halo:nr_halo,ne_wght:ne_wght,np_wght:np_wght,nphi_wght:nphi_wght,$
          emax_wght:emax_wght,ichan_wght:ichan_wght,dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
          wavel_end_wght:wavel_end_wght,calc_npa:calc_npa,calc_spec:calc_spec,calc_birth:calc_birth,calc_fida_wght:calc_fida_wght,$
          calc_npa_wght:calc_npa_wght,calc_brems:calc_brems,load_neutrals:load_neutrals,load_fbm:load_fbm,interactive:interactive}
  
END
