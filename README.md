# FIDASIM
FIDASIM is a code that models the signal that is produced by charge-exchange reactions between fast-ions and injected neutral beams in tokamak plasmas. 
It was originally developed in IDL at UC Irvine [1] and converted to Fortran 90 by Benedikt Geiger (AUGD) [2]. The model is described in Comm. Comp. Phys. 10 (2011) 716. The halo algorithm was improved by Geiger; a weight function calculation was also added.

[1] Heidbrink, W. W., et al. "A code that simulates fast-ion D-alpha and neutral particle measurements." Comm. Comp. Phys. 10 (2011) 716.

[2] Geiger, Benedikt. "Fast-ion transport studies using FIDA spectroscopy at the ASDEX Upgrade tokamak." Diss. lmu, 2013. APA	


***

# How to Install 
## 1. Install dependencies
FIDASIM reads and writes netCDF files. This requires netCDF-4.1.3. You can download the library from [here](http://www.unidata.ucar.edu/downloads/netcdf/netcdf-4_1_3/index.jsp)

Note: By default netCDF will build using the GNU fortran compiler, gfortran. If you plan to use the Intel Fortran Compiler you must also build the 
netCDF library using it. Instructions on how to do this can be found [here](http://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-netcdf-with-the-intel-compilers)
 
## 2. Retrieve FIDASIM source code from GitHub
Clone the git repository from GitHub and change to the source directory: 

    git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
    cd FIDASIM 

## 3. Compile 
FIDASIM will not compile out of the box. You will first need to set the following environmental variables to point to the netCDF lib and include directories. I recommend having this automatically be done at startup.

For tsch shell:

    setenv NETCDF_INCLUDE /path/to/netcdf/install/include
    setenv NETCDF_LIB /path/to/netcdf/install/lib
    setenv LD_LIBRARY_PATH "/path/to/netcdf/install/lib":{$LD_LIBRARY_PATH}
    setenv LD_LIBRARY_PATH "/path/to/netcdf/install/include":{$LD_LIBRARY_PATH}

For bash shell:

    export NETCDF_INCLUDE=/path/to/netcdf/install/include
    export NETCDF_LIB=/path/to/netcdf/install/lib
    export LD_LIBRARY_PATH = /path/to/netcdf/install/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH = /path/to/netcdf/install/include:$LD_LIBRARY_PATH

The Intel Fortran compiler is recommended. You can download the non-commercial version from [here](http://software.intel.com/en-us/non-commercial-software-development)

The last step is the run make in the source directory

    make

## 4. Device Specific Installation Instructions
These installation instructions are unique to each machine. For instructions on how to get FIDASIM to work with a particular machine see the section titled "How do make FIDASIM work for your device"
### DIII-D
FIDASIM currently does not run on the venus cluster since it does not have the required libraries.

### NSTX-U
* To access the git repository and to use the Intel compiler add the following lines to your .login file before the pathscale module is loaded. 

```

module load git/1.8.0.2
module load intel
#GIT AND INTEL MODULES MUST BE LOADED BEFORE PATHSCALE MODULE
module load pathscale

```
* If you run FIDASIM on portal you will get an angry email. Make sure to schedule the job using the "use" command.
* Use the following link to clone the repository http://github.com/D3DEnergetic/FIDASIM.git
* The netCDF directories lib and include, are located at /usr/pppl/intel/11-pkgs/netcdf-4.1.3/

***

#How to run
## 1. Create an input file/procedure.
Prefida, the FIDASIM preprocessing routine, calls an input procedure that returns a structure that contains the input parameters. An example for DIII-D is shown below.
Note: There may be some small differences depending on your device.

```
;;This input file is a procedure so name this file accordingly
PRO input_template,inputs                  ;; Name of this file without .pro

;;-----------------------------------------------------
;;				PREFIDA INPUT FILE
;;-----------------------------------------------------
shot=146088L                               ;; Shot Number
time=1.385                                 ;; Time 
runid='146088H05'                          ;; runid of FIDASIM
device='D3D'                               ;; D3D,NSTX,AUGD,MAST
install_dir='/path/to/install/FIDASIM/'	   ;; Location of fidasim code and executable
result_dir='/path/to/result/directory/'    ;; Location where results will be stored /RESULTS/runid will be made
profile_dir='/path/to/profile/directory/'  ;; Location of profile save files. EX: profile_dir+'shot/'+'dne142353.00505'

;;----------------------------------------------------
;; Fast-ion distribution function from transp
;;----------------------------------------------------
cdf_file='/path/to/the/transp/fast/ion/distribution.cdf'  
                      ;; CDF file from transp with the distribution funciton
emin=0.               ;; Minimum energy used from the distribution function
emax=100.             ;; Maximum energy used from the distribution function
pmin=-1.              ;; Minimum pitch used from the distribution function
pmax=1.	              ;; Maximum pitch used from the distribution function

;;-----------------------------------------------------
;; Beam/FIDA/EQUILIBRIUM Selection
;;-----------------------------------------------------
isource=5             ;; Beam source index (FIDASIM only simulates on NBI source)
einj=0.               ;; [keV] If 0, get data from MDS+
pinj=0.               ;; [MW] If 0, get data from MDS+

diag='OBLIQUE'	      ;; Name of the FIDA diag
equil='EFIT01'        ;; Name of equilibrium. Ex. for D3D EFIT02

;;-----------------------------------------------------
;; Discharge Parameters
;;-----------------------------------------------------
btipsign=-1.d0	      ;; Bt and Ip are in the opposite direction   
ab=2.01410178d0       ;; Atomic mass of beam [u]
ai=2.01410178d0       ;; Atomic mass of hydrogenic plasma ions [u]
impurity_charge=6     ;; 5: BORON, 6: carbon, 7: Nitrogen

;;-----------------------------------------------------
;; Wavelength Grid
;;-----------------------------------------------------
lambdamin=6470.d0                               ;; Minimum wavelength of wavelength grid[A] 
lambdamax=6670.d0                               ;; Maximum wavelength of wavelength grid[A] 
nlambda=2000L                                   ;; Number of wavelengths
dlambda= (lambdamax-lambdamin)/double(nlambda)	;; Wavelength seperation

;;---------------------------------------------------
;; Define FIDASIM grid in machine coordinates(x,y,z)
;;---------------------------------------------------
nx=40                 ;; Number of cells in x direction
ny=60                 ;; Number of cells in y direction
nz=50                 ;; Number of cells in z direction
xdim1=-170.           ;; Minimum x value
xdim2=-70.            ;; Maximum x value
ydim1=-195.           ;; Minimum y value
ydim2=-80.            ;; Maximum y value
zdim1=-70.            ;; Minimum z value
zdim2=70.             ;; Maximum z value

origin=[0.,0.,0.]     ;; If using different a coordinate system, this is the origin 
                      ;; in machine coordinates of the new system

alpha=0.0             ;; Rotation angle in radians from +x about z axis that transforms machine
                      ;; coordinates to the new system. 
beta=0.0              ;; Rotation about +y axis

;;--------------------------------------------------
;; Define number of Monte Carlo particles
;;--------------------------------------------------
nr_fast=5000000       ;; FIDA
nr_ndmc=50000 	      ;; Beam emission
nr_halo=500000        ;; Halo contribution

;;--------------------------------------------------
;; Calculation of the weight function
;;--------------------------------------------------
nr_wght=50            ;; Number of Pitches, energyies and gyro angles 
emax_wght=125.        ;; Maximum energy (keV)
ichan_wght=-1         ;; -1 for all channels, otherwise a given channel index
dwav_wght=.2          ;; Wavelength interval
wavel_start_wght=651. ;; Minimum wavelength
wavel_end_wght=663.   ;; Maximum wavelength

;;-------------------------------------------------
;; Simulation switches
;;-------------------------------------------------
calc_npa=[0]          ;; (0 or 1) If 1 do a simulation for NPA
calc_spec=[1]         ;; (0 or 1) If 1 then spectra is calculated
calc_birth=[1]        ;; (0 or 1) If 1 then the birth profile is calculated
f90brems=[0]          ;; (0 or 1) If 0 use the IDL bremstrahlung calculation
calc_fida_wght=[1]    ;; (0 or 1) If 1 then weight functions are calculated
calc_npa_wght=[0]     ;; (0 or 1) If 1 then weight functions are calculated
load_neutrals=[0]     ;; (0 or 1) If 1 then the neutral density is loaded from an existing 
ps=[0]                ;; (0 or 1) If 1 then make hard copy of plots
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
		wavel_end_wght:wavel_end_wght,calc_npa:calc)npa,calc_spec:calc_spec,calc_birth:calc_birth,calc_fida_wght:calc_fida_wght,$
		calc_npa_wght:calc_npa_wght,f90brems:f90brems,load_neutrals:load_neutrals,ps:ps}

END
```
This template can be found in the TEMPLATES/ directory.

Note: This is an IDL procedure so make sure the procedure name matches the name of the file or else it will not load. 

## Run the pre-processing routine
Prefida pulls in the required profiles and geometry and puts them into netCDF file that FIDASIM can read. Using the above input procedure run the following.

    IDL> prefida,'input_template'

This will make an ASCII input.dat file and and RUNID_inputs.cdf file in a RUNID directory in the result directory. It will also copy the input procedure into the same directory.

Note: prefida can take two keywords: plot and save. 

## Run FIDASIM

    /path/to/fidasim/executable/fidasim /path/to/input/directory/<RUNID>

***

# How do make FIDASIM work for your device
FIDASIM is device agnostic. This means that it doesn't care what your machine is so long as it can read in the needed information.
Prefida, the FIDASIM preprocessing routine, is also device agnostic. It can accomplish this by delegating device-specific routines to other programs. 
In other words, prefida is modular. This allows users to use device-specific routines to prepare the data as long as it delivers the results to prefida in the specified format.
This can best be described with an example. 

Say I have a device called ABCD, which stands for "A Beautiful Cylindrical Device". I want FIDASIM to work with it, so in the source directory I make an ABCD directory.

    mkdir ABCD
    
That directory will hold all of the code that is specific to ABCD. From the TEMPLATES/ directory I copy device_routines.pro into the ABCD directory and rename
it abcd_routines.pro

    cp TEMPLATES/device_routines.pro ABCD/abcd_routines.pro
    
Prefida will look for this routine when it runs with the DEVICE input set to "ABCD". Opening it up I see the following

```
;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO templete_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG

	
	;;IN THIS SECTION YOU CAN USE WHATEVER ROUTINES 
	;;YOU WANT SO LONG AS YOU DEFINE THE OUTPUT STRUCTURES
	;;CONTAIN AT LEAST THE FOLLOWING TAGS

	;;	IDL> help,chords 
	;;	** Structure <1d447c48>, 11 tags, length=728, data length=724, refs=1:
	;;	   NCHAN           LONG                11
	;;	   DIAG            STRING    'OBLIQUE'
	;;	   XLOS            DOUBLE    Array[11]
	;;	   YLOS            DOUBLE    Array[11]
	;;	   ZLOS            DOUBLE    Array[11]
	;;	   XLENS           DOUBLE    Array[11]
	;;	   YLENS           DOUBLE    Array[11]
	;;	   ZLENS           DOUBLE    Array[11]
	;;	   SIGMA_PI_RATIO  DOUBLE    Array[11]
	;;	   HEADSIZE        FLOAT     Array[11]
	;;	   OPENING_ANGLE   FLOAT     Array[11]

	;;	IDL> help,equil
	;;	** Structure <1d474638>, 10 tags, length=6636160, data length=6636138, refs=1:
	;;	   RHO_GRID        FLOAT     Array[40, 60, 50]
	;;	   RHO_CHORDS      STRUCT    -> <Anonymous> Array[1]
	;;	   BX              DOUBLE    Array[40, 60, 50]
	;;	   BY              DOUBLE    Array[40, 60, 50]
	;;	   BZ              DOUBLE    Array[40, 60, 50]
	;;	   EX              DOUBLE    Array[40, 60, 50]
	;;	   EY              DOUBLE    Array[40, 60, 50]
	;;	   EZ              DOUBLE    Array[40, 60, 50]
	;;	   ERR             INT              0

	;;	IDL> help,equil.rho_chords
	;;	** Structure <1d48bf08>, 2 tags, length=352008, data length=352004, refs=2:
	;;	   RHOS            DOUBLE    Array[4000, 11] ;;Rho values along lines of sight
	;;	   DS              FLOAT          0.300000   ;;step size along line of sight in [cm]
	;;

	;;	IDL> help,profiles
	;;	** Structure <1d475698>, 7 tags, length=5816, data length=5810, refs=1:
	;;	   RHO             DOUBLE    Array[121]
	;;	   TE              DOUBLE    Array[121]
	;;	   TI              DOUBLE    Array[121]
	;;	   VTOR            DOUBLE    Array[121]
	;;	   DENE            DOUBLE    Array[121]
	;;	   ZEFF            DOUBLE    Array[121]
	;;	   ERR             INT              0

	;;	IDL> help,nbi
	;;	** Structure <1d475af8>, 13 tags, length=168, data length=168, refs=1:
	;;	   EINJ            DOUBLE           80.775734
	;;	   PINJ            DOUBLE           2.4117758
	;;	   FULL            DOUBLE          0.54850105
	;;	   HALF            DOUBLE          0.28972649
	;;	   THIRD           DOUBLE          0.16177245
	;;	   XYZ_SRC         DOUBLE    Array[3]
	;;	   XYZ_POS         DOUBLE    Array[3]
	;;	   BMWIDRA         DOUBLE           6.0000000
	;;	   BMWIDZA         DOUBLE           24.000000
	;;	   DIVY            DOUBLE    Array[3]
	;;	   DIVZ            DOUBLE    Array[3]
	;;	   FOCY            DOUBLE           999999.90
	;;	   FOCZ            DOUBLE           1000.0000

	;;FOR CONVENIENCE HERE ARE THE MINIMUM STRUCTURE DEFINITIONS
	equil={rho_grid:rho_grid,$	   			;;FIDA GRID IN MAGNETIC FLUX COORDINATES (RHO)
		   rho_chords:rho_chords,$			;;STRUCTURE CONTAINING AN ARRAY OF RHO VALUES AND STEP SIZE IN [cm]
		   bx:bx,$					   		;;X MAGNETIC FIELD COMPONENT AT GRID POINTS
		   by:by,$					   		;;Y MAGNETIC FIELD COMPONENT AT GRID POINTS
		   bz:bz,$					   		;;Z MAGNETIC FIELD COMPONENT AT GRID POINTS
		   ex:ex,$							;;X ELECTRIC FIELD COMPONENT AT GRID POINTS
		   ey:ey,$							;;Y ELECTRIC FIELD COMPONENT AT GRID POINTS
		   ez:ez }							;;Z ELECTRIC FIELD COMPONENT AT GRID POINTS

	nbi={einj:einj,$				   		;;BEAM INJECTION ENERGY [keV]
		 pinj:pinj,$				   		;;BEAM INJECTION POWER  [MW]
		 full:full,$				   		;;FULL BEAM FRACTION
		 half:half,$				   		;;HALF BEAM FRACTION
		 third:third,$				   		;;THIRD BEAM FRACTION
		 xyz_src:xyz_src,$			   		;;POSITION OF BEAM SOURCE IN MACHINE COORDINATES [cm]
		 xyz_pos:xyz_pos,$			   		;;BEAM CROSSOVER POINT IN MACHINE COORDINATES [cm]
		 bmwidra:bmwidra,$			   		;;HORIZONTAL BEAM WIDTH [cm]
		 bmwidza:mbwidza,$			   		;;VERTICAL BEAM WIDTH   [cm]
		 focy:focy,$				   		;;HORIZONTAL FOCAL LENGTH [cm]
	     focz:focz,$						;;VERTICAL FOCAL LENGTH [cm]
		 divy:divy,$				   		;;HORIZONTAL BEAM DIVERGENCE [rad]
		 divz:divz }				   		;;VERTICAL BEAM DIVERGENCE [rad]

  	chords={sigma_pi_ratio:sigma_pi_ratio,$	;;RATIO OF SIGMA LINES TO PI LINES
		 nchan:nchan,$				  		;;NUMBER OF CHANNELS
		 xmid:xmid,$						;;X POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 ymid:ymid,$						;;Y POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
         zmid:zmid,$						;;Z POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 xlens:xlens,$						;;X POS. OF LENS [cm]
		 ylens:ylens,$						;;Y POS. OF LENS [cm]
 		 zlens:zlens,$						;;Z POS. OF LENS [cm]
  		 headsize:headsize,$				;;SIZE OF HEAD
		 opening_angle:opening_angle}		;;OPENING ANGLE

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  vtor:vtor,$					;;TORODIAL ANGULAR VELOCITY [rad/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END 
```
As you can see it is all rather self explaintory. All I need to do now to get ABCD to work with prefida is supply the above information.
I can use whatever routines I want so long as the output of abcd_routines is defined as above. I could also include other things in the
structures so long as I have the above. 

Note: All positions in the device routines are in machine coordinates

# Frequently Asked Questions
### I get a segmentation fault when I run the code with multiple cores, but not when I only use one core. 
A segmentation fault happens when a code tries to access memory that it wasn't allocated. When a code runs in parallel each thread is allocated a certain amount of memory. Under certain conditions, i.e. a large grid, a thread runs out of room and tries to access more. This throws a seg fault. To increase this limit run the following:

For tsch shell:

    limit stacksize unlimited
For bash shell:

    ulimit -s unlimited
    
### I have fixed your code. How do I submit a patch?
FIDASIM is open source code. In order to contribute to the project please fork us on GitHub and open a pull request. We will then review your changes and incorporate them into code base.

### The plot keyword doesn't do anything. 
Plotting, like the device routines, are often very specific to the device. Accordingly, plotting is treated the same way as the device routines. 
Prefida will look for a routine called, using the nomenclature from above, abcd_plots.pro. A basic plot routine can be found in the TEMPLATES/ directory.
