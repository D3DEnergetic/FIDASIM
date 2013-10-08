PRO east_input,inputs                                   
;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=34616		;; Shot Number
time=3.5 	   	;; Time (s) 
runid='weight_test'   	;; runid of FIDASIM
device='EAST'										
install_dir='/home/FIDA/FIDASIM/'		   ;; Location of fidasim code and executable
result_dir=install_dir+'RESULTS/EAST/'  	;; Location where results will be stored result_dir/runid directory will be created
transpid='34616B02'             ;; Name of transp run (used for plasma parameters)
profile_dir=install_dir+'EAST/data'                  ;; Directory with transp *.CDF file

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/home/FIDA/FIDASIM/EAST/data/34616B02_fi_1.cdf'  ;; CDF file from transp with the distribution function
emin=0.                    ;; minimum energy (keV) used from the distribution function
emax=100.        	 ;; maximum energy (keV) used from the distribution function
pmin=-1.          ;; minimum pitch used from the distribution function
pmax=1.          ;; maximum pitch used from the distribution function
;;-----------------------------------------------------
;; Beam/diagnostic/equilibrium Selection
;;-----------------------------------------------------
isource=0     		    ;; Beam source index (FIDASIM only simulates one NBI source)
einj=80.                 ;; [keV] If 0, retrieves archived data (not implemnented)
pinj=4.0                 ;; [MW] If 0, retrieves archived data

diag='VERTICAL-B'	          ;; Name of the FIDA or NPA diagnostic

gfile='/home/FIDA/FIDASIM/EAST/data/g034616.03200'   ;; EFIT g eqdsk

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
origin=[0.,0.,0.0]
alpha=0.D 	;; rotation about z axis
beta=0.0
nx=25				;; Number of cells in x direction
ny=30				;; Number of cells in y direction
nz=25				;; Number of cells in z direction
xdim1=150.			;; Minimum x value
xdim2=230.			;; Maximum x value
ydim1=0.			;; Minimum y value
ydim2=90.			;; Maximum y value
zdim1=-40.			;; Minimum z value
zdim2=40.			;; Maximum z value

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
;nr_fast=5000000   		;; FIDA
;nr_ndmc=100000 			;; Beam emission
;nr_halo=1000000	   		;; Halo contribution
nr_fast=50000   		;; FIDA
nr_ndmc=10000 			;; Beam emission
nr_halo=10000	   		;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
nr_wght=30   				;; Number of Pitches, energyies and gyro angles 
emax_wght=90  				;; Maximum energy (keV)
ichan_wght=-1.  				;; -1 for all channels, otherwise a given channel index
dwav_wght=1.   				;; Wavelength interval
wavel_start_wght=650.1  		;; Minimum wavelength
wavel_end_wght=661.1   		;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[0]   				;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[0]   			;; (0 or 1) If 1 then spectra are calculated
sim_fida=[0]    			;; (0 or 1) If 1 then the FIDA spectra are simulated
calc_birth=[0]                          ;; (0 or 1) If 1 then the birth profile is calculated
f90brems=[0]                ;; (0 or 1) If 0 use the IDL bremsstrahlung calculation
calc_fida_wght=[1]  		;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[0]  			;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]   	;; (0 or 1) If 1 then the neutral density is loaded from an existing 
ps=[0]                  ;; (0 or 1) If 1 then make hard copy of plots

;;------------------------------------------------
;; DO NOT MODIFY THIS PART
;;------------------------------------------------

inputs={shot:shot,time:time,runid:runid,device:strupcase(device),install_dir:install_dir, $
        result_dir:result_dir,transpid:transpid,profile_dir:profile_dir, $
       cdf_file:cdf_file,emin:emin,emax:emax,pmin:pmin,pmax:pmax, $
       isource:isource,einj:einj,pinj:pinj,diag:diag,gfile:gfile,$
       btipsign:btipsign,ab:ab,ai:ai,impurity_charge:impurity_charge,$
       lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,origin:origin,alpha:alpha,beta:beta,$
       nx:nx,ny:ny,nz:nz,xdim1:xdim1,xdim2:xdim2,ydim1:ydim1,ydim2:ydim2,zdim1:zdim1,zdim2:zdim2,$
       nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo,nr_wght:nr_wght,$
       emax_wght:emax_wght,ichan_wght:ichan_wght,dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
	   wavel_end_wght:wavel_end_wght,calc_npa:calc_npa,calc_spec:calc_spec, $
       sim_fida:sim_fida, calc_birth:calc_birth, $
       f90brems:f90brems,calc_fida_wght:calc_fida_wght,calc_npa_wght:calc_npa_wght, $
       load_neutrals:load_neutrals,ps:ps}

END
