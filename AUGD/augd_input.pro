;;This input file is a procedure so name this file accordingly
PRO augd_input,inputs                                   ;; Name of this file without .pro
user=+GETENV('USER')
;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=28746L 											;; Shot Number
time=4.421  											;; Time 
runid='28746A01'	   									;; runid of FIDASIM
device='AUGD'											;; D3D,NSTX,AUGD,MAST
install_dir='/u/'+user+'/FIDASIM/'						;; Location of fidasim code and executable
result_dir='/u/'+user+'/FIDASIM/RESULTS/AUGD/'  		;; Location where results will be stored /RESULTS/runid will be made
profile_dir='/u/'+user+'/GAPROFILES/'					;; Location of profile save files. EX: profile_dir+'shot/'+'dne142353.00505'

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/u/'+user+'/FIDASIM/TEST/28746A01_fi_1.cdf'        ;; CDF file from transp with the distribution funciton
emin=0.    													 ;; Minimum energy used from the distribution function
emax=100.  													 ;; Maximum energy used from the distribution function
pmin=-1.													 ;; Minimum pitch used from the distribution function
pmax=1.														 ;; Maximum pitch used from the distribution function

;;-----------------------------------------------------
;; Beam/FIDA/EQUILIBRIUM Selection
;;-----------------------------------------------------
isource=2     											;; Beam source index (FIDASIM only simulates on NBI source)
einj=0.                 								;; [keV] If 0, get data from MDS+
pinj=0.                 								;; [MW] If 0, get data from MDS+

diag='CFR'												;; Name of the FIDA diag
equil='EQH'												;; Name of equilibrium. Ex. for D3D EFIT02

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=-1.d0											;; Bt and Ip are in the opposite direction   
ab=2.01410178d0             							;; Atomic mass of beam [u]
ai=2.01410178d0             							;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=5           							;; 5: BORON, 6: carbon, 7: Nitrogen

;;-----------------------------------------------------
;; Wavelength Grid
;;-----------------------------------------------------
lambdamin=6470.d0           							;; Minimum wavelength of wavelength grid[A] 
lambdamax=6670.d0           							;; Maximum wavelength of wavelength grid[A] 
nlambda=2000L               							;; Number of wavelengths
dlambda= (lambdamax-lambdamin)/double(nlambda)			;; Wavelength seperation

;;---------------------------------------------------
;; Define FIDASIM grid in machine coordinates(x,y,z)
;;---------------------------------------------------
nx=40               ;; Number of cells in x direction
ny=60               ;; Number of cells in y direction
nz=50               ;; Number of cells in z direction
xdim1=100.          ;; Minimum x value
xdim2=220.          ;; Maximum x value
ydim1=-150.         ;; Minimum y value
ydim2=40.           ;; Maximum y value
zdim1=-75.          ;; Minimum z value
zdim2=75.           ;; Maximum z value

origin=[0.,0.,0.]	;; If using different a coordinate system, this is the origin 
					;; in machine coordinates of the new system

alpha=0.0		    ;; Rotation angle in radians from +x about z axis that transforms machine
					;; coordinates to the new system. Rotations about x or y not supported
beta=0.0
;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fast=500000   										;; FIDA
nr_ndmc=5000 											;; Beam emission
nr_halo=50000   										;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
nr_wght=50   											;; Number of Pitches, energyies and gyro angles 
emax_wght=100.  										;; Maximum energy (keV)
ichan_wght=-1  											;; -1 for all channels, otherwise a given channel index
dwav_wght=.2   											;; Wavelength interval
wavel_start_wght=651.  									;; Minimum wavelength
wavel_end_wght=663.   									;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[0]   											;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[1]   										;; (0 or 1) If 1 then spectra is calculated
calc_birth=[1]    										;; (0 or 1) If 1 then the birth profile is calculated
f90brems=[1]                							;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_fida_wght=[1]  									;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[0]  										;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]   									;; (0 or 1) If 1 then the neutral density is loaded from an existing 
														;; neutrals.bin file located in runid directory
ps=[0]													;; (0 or 1) If 1 then make hard copy of plots
;;------------------------------------------------
;; DO NOT MODIFY THIS PART
;;------------------------------------------------

inputs={shot:shot,time:time,runid:runid,device:strupcase(device),install_dir:install_dir,result_dir:result_dir,$
	    cdf_file:cdf_file,profile_dir:profile_dir,emin:emin,emax:emax,pmin:pmin,pmax:pmax,isource:isource,diag:diag,$
	    einj:einj,pinj:pinj,equil:equil,btipsign:btipsign,ab:ab,ai:ai,impurity_charge:impurity_charge,$
	    lambdamin:lambdamin,lambdamax:lambdamax,nlambda:nlambda,dlambda:dlambda,$
	    nx:nx,ny:ny,nz:nz,xdim1:xdim1,xdim2:xdim2,ydim1:ydim1,ydim2:ydim2,zdim1:zdim1,zdim2:zdim2,$
		origin:origin,alpha:alpha,beta:beta,nr_fast:nr_fast,nr_ndmc:nr_ndmc,nr_halo:nr_halo,nr_wght:nr_wght,$
        emax_wght:emax_wght,ichan_wght:ichan_wght,dwav_wght:dwav_wght,wavel_start_wght:wavel_start_wght,$
		wavel_end_wght:wavel_end_wght,calc_npa:calc_npa,calc_spec:calc_spec,calc_birth:calc_birth,calc_fida_wght:calc_fida_wght,$
		calc_npa_wght:calc_npa_wght,f90brems:f90brems,load_neutrals:load_neutrals,ps:ps}

END
