# FIDASIM
FIDASIM is a code that models the signal that is produced by charge-exchange reactions between fast-ions and injected neutral beams in tokamak plasmas. 

***

# How to Install 
## 1. Install dependencies
FIDASIM reads and writes netCDF files. In order to use this code you must download the required fortran libraries. You can download the library from [here](https://github.com/Unidata/netcdf-fortran/releases)

## 2. Retrieve FIDASIM source code from GitHub
Clone the git repository from GitHub and change to the source directory: 

    git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
	cd FIDASIM 

## 3. Compile 
FIDASIM will not compile out of the box. You will first need to set the following environmental variables to point to the netCDF lib and include directories. I recommend having this automatically be done at startup.
For tsch:
    setenv NETCDF_INCLUDE /path/to/netcdf/install/include
    setenv NETCDF_LIB /path/to/netcdf/install/lib
For bash:
    export NETCDF_INCLUDE=/path/to/netcdf/install/include
    export NETCDF_LIB=/path/to/netcdf/install/lib

The Intel Fortran compiler is recommended. You can download the non-commercial version from [here](http://software.intel.com/en-us/non-commercial-software-development)

The last step is the run make in the source directory
    make

## 4. Device Specific Notes
### DIII-D
FIDASIM currently does not run on the venus cluster since it does not have the required libraries.

### NSTX-U
* To access the git repository and to use the Intel compiler add the following lines to your .login file before the pathscale module is loaded. 
	module load git
    module load intel
	#GIT AND INTEL MODULES MUST BE LOADED BEFORE PATHSCALE MODULE
	module load pathscale
* If you run FIDASIM on portal you will get an angry email. Make sure to schedule the job using the "use" command.

***

#How to run
## 1. Create an input file/procedure.
prefida, the FIDASIM preprocessing routine, calls an input procedure that returns a structure that contains the input parameters. An example is shown below.

''''
;;This input file is a procedure so name this file accordingly
PRO input_template,inputs                                   ;; Name of this file without .pro

;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=146088L											;; Shot Number
time=1.385  											;; Time 
runid='146088H05'   									;; runid of FIDASIM
device='D3D'											;; D3D,NSTX,AUGD,MAST
install_dir='/path/to/install/FIDASIM/'					;; Location of fidasim code and executable
result_dir='/path/to/result/directory/'  				;; Location where results will be stored /RESULTS/runid will be made
profile_dir='/path/to/profile/directory/'				;; Location of profile save files. EX: profile_dir+'shot/'+'dne142353.00505'

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/path/to/the/transp/fast/ion/distribution.cdf'     ;; CDF file from transp with the distribution funciton
emin=0.    													 ;; Minimum energy used from the distribution function
emax=100.  													 ;; Maximum energy used from the distribution function
pmin=-1.													 ;; Minimum pitch used from the distribution function
pmax=1.														 ;; Maximum pitch used from the distribution function

;;-----------------------------------------------------
;; Beam/FIDA/EQUILIBRIUM Selection
;;-----------------------------------------------------
isource=5     											;; Beam source index (FIDASIM only simulates on NBI source)
einj=0.                 								;; [keV] If 0, get data from MDS+
pinj=0.                 								;; [MW] If 0, get data from MDS+

diag='OBLIQUE'											;; Name of the FIDA diag
equil='EFIT01'											;; Name of equilibrium. Ex. for D3D EFIT02

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=-1.d0											;; Bt and Ip are in the opposite direction   
ab=2.01410178d0             							;; Atomic mass of beam [u]
ai=2.01410178d0             							;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=6           							;; 5: BORON, 6: carbon, 7: Nitrogen

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
xdim1=-170.         ;; Minimum x value
xdim2=-70.          ;; Maximum x value
ydim1=-195.         ;; Minimum y value
ydim2=-80.          ;; Maximum y value
zdim1=-70.          ;; Minimum z value
zdim2=70.           ;; Maximum z value

origin=[0.,0.,0.]	;; If using different a coordinate system, this is the origin 
					;; in machine coordinates of the new system

alpha=0.0		    ;; Rotation angle in radians from +x about z axis that transforms machine
					;; coordinates to the new system. 
beta=0.0			;; Rotation about +y axis

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fast=5000000   										;; FIDA
nr_ndmc=50000 											;; Beam emission
nr_halo=500000   										;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
nr_wght=50   											;; Number of Pitches, energyies and gyro angles 
emax_wght=125.  										;; Maximum energy (keV)
ichan_wght=-1  											;; -1 for all channels, otherwise a given channel index
dwav_wght=.2   											;; Wavelength interval
wavel_start_wght=651.  									;; Minimum wavelength
wavel_end_wght=663.   									;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
npa=[0]   												;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[1]   										;; (0 or 1) If 1 then spectra is calculated
sim_fida=[1]											;; (0 or 1) if 1 then the FIDA spectra is calculated
calc_birth=[1]    										;; (0 or 1) If 1 then the birth profile is calculated
f90brems=[0]                							;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_wght=[1]  											;; (0 or 1) If 1 then weight functions are calculated
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
		wavel_end_wght:wavel_end_wght,npa:npa,calc_spec:calc_spec,sim_fida:sim_fida,calc_birth:calc_birth,calc_wght:calc_wght,$
		f90brems:f90brems,load_neutrals:load_neutrals,ps:ps}

END
''''
This template can be found in the TEMPLATES/ directory.

Note: This is an IDL procedure so make sure the procedure name matches the name of the file or else it will not load. 

## Run the pre-processing routine
Prefida pulls in the required profiles and geometry and puts them into netCDF file that FIDASIM can read. Using the above input procedure run the following.
    IDL> prefida,'input_template'
This will make an ASCII input.dat file and and <RUNID>_inputs.cdf file in a <RUNID> directory in the result directory. It will also copy the input procedure into the same directory.

Note: prefida can take two keywords: plot and save. 

## Run FIDASIM
    /path/to/fidasim/executable/fidasim /path/to/input/directory/<RUNID>


# References

Heidbrink, W. W., et al. "A code that simulates fast-ion D-alpha and neutral particle measurements." Comm. Comp. Phys. 10 (2011) 716.

Geiger, Benedikt. "Fast-ion transport studies using FIDA spectroscopy at the ASDEX Upgrade tokamak." Diss. lmu, 2013. APA	

