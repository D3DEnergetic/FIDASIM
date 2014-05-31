PRO d3d_input,inputs                                   
user=+GETENV('USER')
;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=142114											;; Shot Number
time=1.090 		        							;; Time (s) 
runid='142114B08'   								;; runid of FIDASIM
device='D3D'										;; D3D,NSTX,AUGD,MAST
result_dir='/u/'+user+'/FIDASIM/RESULTS/D3D/'  	;; Location where results will be stored result_dir/runid directory will be created
profile_dir='/u/'+user+'/GAPROFILES/'              	;; Location of profile save files. EX: profile_dir+'shot/'+'dne142353.00505'

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/u/'+user+'/GAPROFILES/142114/142114B08_fi_1.cdf'  ;; CDF file from transp with the distribution funciton
emin=0.                                                      ;; minimum energy used from the distribution function
emax=100.        										     ;; maximum energy used from the distribution function
pmin=-1.
pmax=1.

;;-----------------------------------------------------
;; Beam/diagnostic/equilibrium Selection
;;-----------------------------------------------------
isource=6     		    ;; Beam source index (FIDASIM only simulates one NBI source)
einj=0.                 ;; [keV] If 0, get data from MDS+
pinj=0.                 ;; [MW] If 0, get data from MDS+

diag='VERTICAL'	;; Name of the diagnostic

equil='EFIT01'			;; MDS+ equilibrium runid

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=-1.d0		;; Use -1 when Bt and Ip are in the opposite direction   
ab=2.01410178d0     ;; Atomic mass of beam [u]
ai=2.01410178d0     ;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=6 	;; 5: BORON, 6: carbon, 7: Nitrogen

;;-----------------------------------------------------
;; Wavelength Grid
;;-----------------------------------------------------
lambdamin=6470.d0       						;; Minimum wavelength of wavelength grid[A] 
lambdamax=6670.d0       						;; Maximum wavelength of wavelength grid[A] 
nlambda=2000L           						;; Number of wavelengths
dlambda= (lambdamax-lambdamin)/double(nlambda)	;; Wavelength seperation

;;---------------------------------------------------
;; Define FIDASIM grid in machine coordinates(x,y,z)
;;---------------------------------------------------
origin=[-180.16373,600.82700,0.0]
alpha=-4.6020171774D 	;; rotation about z axis
beta=0.0
nx=40				;; Number of cells in x direction
ny=60				;; Number of cells in y direction
nz=50				;; Number of cells in z direction
xmin=-520.			;; Minimum x value
xmax=-420.			;; Maximum x value
ymin=-40.			;; Minimum y value
ymax=40.			;; Maximum y value
zmin=-40.			;; Minimum z value
zmax=40.			;; Maximum z value

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fida=5000   			;; FIDA
nr_ndmc=1000 			;; Beam emission
nr_halo=5000	   		;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
ne_wght=50                  ;; Number of Energies
np_wght=50                  ;; Number of Pitches
nphi_wght=50   				;; Number of Gyro-angles 
emax_wght=100  				;; Maximum energy (keV)
ichan_wght=-1  				;; -1 for all channels, otherwise a given channel index
dwav_wght=1.   				;; Wavelength interval
wavel_start_wght=651.  		;; Minimum wavelength
wavel_end_wght=663.   		;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[0]   				;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[0]    		    ;; (0 or 1) If 1 then spectra is calculated
calc_birth=[0]              ;; (0 or 1) If then the birth profile is calculated
f90brems=[0]                ;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_fida_wght=[0]  		;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[0]  			;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]   		;; (0 or 1) If 1 then the neutral density is loaded from an existing 

;;------------------------------------------------
;; DO NOT MODIFY THIS PART
;;------------------------------------------------
install_dir=+GETENV('FIDASIM_DIR')
inputs={shot:shot,time:time,runid:runid,device:strupcase(device),install_dir:install_dir,result_dir:result_dir,$
       cdf_file:cdf_file,profile_dir:profile_dir,emin:emin,emax:emax,pmin:pmin,pmax:pmax, $
       isource:isource,einj:einj,pinj:pinj,diag:diag,equil:equil,$
       btipsign:btipsign,ab:ab,ai:ai,impurity_charge:impurity_charge,$
       lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,origin:origin,alpha:alpha,beta:beta,$
       nx:nx,ny:ny,nz:nz,xmin:xmin,xmax:xmax,ymin:ymin,ymax:ymax,zmin:zmin,zmax:zmax,$
       nr_fida:nr_fida,nr_ndmc:nr_ndmc,nr_halo:nr_halo,ne_wght:ne_wght,np_wght:np_wght,nphi_wght:nphi_wght,$
       emax_wght:emax_wght,ichan_wght:ichan_wght,dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
	   wavel_end_wght:wavel_end_wght,calc_npa:calc_npa,calc_spec:calc_spec,calc_birth:calc_birth, $
       f90brems:f90brems,calc_fida_wght:calc_fida_wght,calc_npa_wght:calc_npa_wght,load_neutrals:load_neutrals}

END
