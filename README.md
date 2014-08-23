# FIDASIM
FIDASIM is a code that models the signal that is produced by charge-exchange reactions between fast-ions and injected neutral beams in tokamak plasmas. 
It was originally developed in IDL at UC Irvine [1] and converted to Fortran 90 by Benedikt Geiger (AUGD) [2]. The model is described in Comm. Comp. Phys. 10 (2011) 716. The halo algorithm was improved by Geiger; a weight function calculation was also added.

[1] Heidbrink, W. W., et al. "A code that simulates fast-ion D-alpha and neutral particle measurements." Comm. Comp. Phys. 10 (2011) 716.

[2] Geiger, Benedikt. "Fast-ion transport studies using FIDA spectroscopy at the ASDEX Upgrade tokamak." Diss. lmu, 2013. APA	


***

# How to Install 
## 1. Install dependencies
FIDASIM reads and writes netCDF files. This requires netCDF-4.1.3 or earlier. You can download the library from [here](http://www.unidata.ucar.edu/downloads/netcdf/netcdf-4_1_3/index.jsp)

Note: By default netCDF will build using the GNU Fortran compiler, gfortran. If you plan to use the Intel Fortran compiler you must also build the 
netCDF library using it. Instructions on how to do this can be found [here](http://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-netcdf-with-the-intel-compilers).
Also, netCDF has the option of using HDF5 data format. This, naturally, requires the HDF5 libraries. If you do not have access to the HDF5 libraries netCDF can be built without it.

## 2. Retrieve FIDASIM source code from GitHub
Clone the git repository from GitHub and change to the source directory: 

    git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
    cd FIDASIM 

## 3. Compile 
FIDASIM will not compile out of the box. You will first need to set the following environmental variables to point to the netCDF and install directories. I recommend having this automatically be done at startup.

For tsch shell:

    setenv FIDASIM_DIR /path/to/fidasim/install/    #don't forget the last slash
    setenv FIDASIM_COMPILER gfortran #use 'ifort' for Intel compiler
    setenv NETCDF_INCLUDE /path/to/netcdf/install/include
    setenv NETCDF_LIB /path/to/netcdf/install/lib
    setenv LD_LIBRARY_PATH "/path/to/netcdf/install/lib":{$LD_LIBRARY_PATH}
    setenv LD_LIBRARY_PATH "/path/to/netcdf/install/include":{$LD_LIBRARY_PATH}
    setenv PATH {$FIDASIM_DIR}LIB:{$PATH}
    
For bash shell:

    export FIDASIM_DIR=/path/to/fidasim/install/    #don't forget the last slash
    export FIDASIM_COMPILER=gfortran #use 'ifort' for Intel compiler
    export NETCDF_INCLUDE=/path/to/netcdf/install/include
    export NETCDF_LIB=/path/to/netcdf/install/lib
    export LD_LIBRARY_PATH=/path/to/netcdf/install/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/path/to/netcdf/install/include:$LD_LIBRARY_PATH
    export PATH=$FIDASIM_DIR/LIB:$PATH

The GNU Fortran compiler, gfortran (version >= 4.7.3), is recommended. You can also use Intel Fortran compiler (version >= 11.0). The non-commercial version of the Intel compiler is available [here](http://software.intel.com/en-us/non-commercial-software-development)

The last step is the run make in the source directory

    make

## 4. Run a test case
From the install directory run 

    fidasim TEST/D3D/146088H06_inputs.dat

## 5. Device Specific Installation Instructions
These installation instructions are unique to each machine. For instructions on how to get FIDASIM to work with a particular machine see the section titled "How do make FIDASIM work for your device"
### DIII-D
Cerview routines are needed to run prefida. Add the commands located in ```D3D/d3d_startup.pro``` to your IDL startup file or start IDL as follows
```bash
venus ~ $ idl D3D/d3d_startup.pro
```
 
FIDASIM currently does not run on the venus cluster since it does not have the required libraries.

### NSTX-U
* To access the git repository and to use the Intel compiler add the following lines to your .login file before the pathscale module is loaded. 

```

module load git/1.8.0.2
module load intel
#GIT AND INTEL MODULES MUST BE LOADED BEFORE PATHSCALE MODULE
module load pathscale
module load nstx/python-2.7
module load python/scipy
```
* If you run FIDASIM on portal you will get an angry email. Make sure to schedule the job using the "use" command.
* Use the following link to clone the repository http://github.com/D3DEnergetic/FIDASIM.git
* The netCDF directories lib and include, are located at /usr/pppl/intel/11-pkgs/netcdf-4.1.3/

***

#How to run
## 1. Create an input file
Prefida, the FIDASIM preprocessing routine, reads in a JSON file that contains the input parameters. An example for DIII-D is shown below.
Note: There may be some small differences depending on your device.
Note: JSON input file only runs on IDL v8.2 and above. If you don't have a compatible version see `TEMPLATES/input_template.pro` 
```C
//-----------------------------------------------------
//                PREFIDA INPUT FILE
//-----------------------------------------------------
{
"shot":146088,          // Shot Number
"time":1.385,           // Time 
"runid":"146088H05",    // runid of FIDASIM
"device":"D3D",         // D3D",NSTX",AUGD",MAST
"result_dir":"/u/stagnerl/FIDASIM/RESULTS/D3D/",           
                        // Location where results will be stored
"profile_dir":"/u/heidbrin/OANB/AUG/",                     
                        // Location of profile files

//-----------------------------------------------------
// Fast-ion distribution function from transp
//-----------------------------------------------------
"cdf_file":"/e/alfven/FIDAsim/D3D/146088/146088H02_fi_9.cdf",    
                        // CDF file from transp with the distribution funciton
"emin":0.0,             // Minimum energy used from the distribution function
"emax":100.0,           // Maximum energy used from the distribution function
"pmin":-1.0,            // Minimum pitch used from the distribution function
"pmax":1.0,             // Maximum pitch used from the distribution function

//-----------------------------------------------------
// Beam/FIDA/EQUILIBRIUM Selection
//-----------------------------------------------------
"isource":5,            // Beam source index (FIDASIM only simulates one NBI source)
"einj":0.,              // [keV] If 0, get data from MDS+
"pinj":0.,              // [MW] If 0, get data from MDS+
"diag":"OBLIQUE",       // Name of the FIDA diag
"equil":"EFIT01",       // Name of equilibrium. Ex. for D3D EFIT02

//-----------------------------------------------------
// Discharge Parameters
//-----------------------------------------------------
"btipsign":-1.0,        // Bt and Ip are in the opposite direction   
"ab":2.01410178,        // Atomic mass of beam [u]
"ai":2.01410178,        // Atomic mass of hydrogenic plasma ions [u]
"impurity_charge":6,    // 5: BORON, 6: carbon, 7: Nitrogen

//-----------------------------------------------------
// Wavelength Grid
//-----------------------------------------------------
"lambdamin":647.0,      // Minimum wavelength of wavelength grid [nm] 
"lambdamax":667.0,      // Maximum wavelength of wavelength grid [nm] 
"nlambda":2000,         // Number of wavelengths
"dlambda":0.01,         // Wavelength seperation

//---------------------------------------------------
// Define FIDASIM grid in machine coordinates(x,y,z)
//---------------------------------------------------
"nx":40,                // Number of cells in x direction
"ny":60,                // Number of cells in y direction
"nz":50,                // Number of cells in z direction
"xmin":-170.0,         // Minimum x value
"xmax":-70.0,          // Maximum x value
"ymin":-195.0,         // Minimum y value
"ymax":-80.0,          // Maximum y value
"zmin":-70.0,          // Minimum z value
"zmax":70.0,           // Maximum z value

"origin":[0,0,0],       // If using different a coordinate system, this is the origin 
                        // in machine coordinates of the new system

"alpha":0.0,            // Rotation angle in radians from +x about z axis that transforms machine
                        // coordinates to the new system. 
"beta":0.0,             // Rotation about +y axis

//--------------------------------------------------
// Define number of Monte Carlo particles
//--------------------------------------------------
"nr_fast":5000000,      // FIDA
"nr_nbi":50000,         // Beam emission
"nr_halo":500000,       // Halo contribution

//--------------------------------------------------
// Calculation of the weight function
//--------------------------------------------------
"ne_wght":50,               // Number of Energies 
"np_wght":50,               // Number of Pitches 
"nphi_wght":50,             // Number of Gyro-angles 
"emax_wght":125.0,          // Maximum energy (keV)
"ichan_wght":-1,            // -1 for all channels", otherwise a given channel index
"dwav_wght":0.2,            // Wavelength interval
"wavel_start_wght":651.0,   // Minimum wavelength
"wavel_end_wght":663.0,     // Maximum wavelength

//-------------------------------------------------
// Simulation switches
//-------------------------------------------------
"calc_npa":0,           // (0 or 1) If 1 do a simulation for NPA
"calc_spec":1,          // (0 or 1) If 1 then spectra is calculated
"calc_birth":1,         // (0 or 1) If 1 then the birth profile is calculated
"calc_brems":0,         // (0 or 1) If 0 use the IDL bremstrahlung calculation
"calc_fida_wght":1,     // (0 or 1) If 1 then fida weight functions are calculated
"calc_npa_wght":0,      // (0 or 1) If 1 then npa weight functions are calculated
"load_neutrals":0,      // (0 or 1) If 1 then the neutral density is loaded from an existing run 
"load_fbm":1            // (0 or 1) If 1 then the fbm is loaded (calc_spec/npa overwrites)

//------------------------------------------------
// Extra Variables
//------------------------------------------------

}
```
This template can be found in the TEMPLATES/ directory.

Note: The input file is not proper JSON since it has C style comments. The comments are stripped out of the file before being parsed.

## Run the pre-processing routine
Prefida pulls in the required profiles and geometry and puts them into netCDF file that FIDASIM can read. Using the above input file run the following.

    IDL> prefida,'input_template.json'

This will make an FORTRAN namelist file and a netCDF inputs file in the result directory. It will also copy the input file into the same directory.

Note: prefida can take two keywords: plot and save. 

## Run FIDASIM

    /path/to/fidasim/executable/fidasim /path/to/input/directory/<RUNID>_inputs.dat

## Read the output files
FIDASIM can create the following output netCDF files depending on the simulation switches
```
<RUNID>_neutrals.cdf: This file contains the neutral beam and halo density.
<RUNID>_birth.cdf: This file contains the birth positions of the fast ions.
<RUNID>_spectra.cdf: This file contains the simulated spectra.
<RUNID>_npa.cdf: This file contains the simulated NPA spectrum.
<RUNID>_fida_weights.cdf: This file contains the calculated FIDA weight functions.
<RUNID>_npa_weights.cdf: This file contains the calculated NPA flux and weight functions.
```
The netCDF files can be read by using the IDL routine ```read_ncdf.pro``` or ```load_results.pro``` located in the ```LIB``` directory
***

# How do make FIDASIM work for your device
FIDASIM is device agnostic. This means that it doesn't care what your machine is so long as it can read in the needed information.
Prefida, the FIDASIM preprocessing routine, is also device agnostic. It can accomplish this by delegating device-specific routines to other programs. 
In other words, prefida is modular. This allows users to use device-specific routines to prepare the data as long as it delivers the results to prefida in the specified format.
This can best be described with an example. 

Say I have a device called ABCD, which stands for "A Big Cylindrical Device". I want FIDASIM to work with it, so in the source directory I make an ABCD directory.

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
	;;	   CHAN_ID         LONG      Array[11]
	;;	   XLOS            DOUBLE    Array[11]
	;;	   YLOS            DOUBLE    Array[11]
	;;	   ZLOS            DOUBLE    Array[11]
	;;	   XLENS           DOUBLE    Array[11]
	;;	   YLENS           DOUBLE    Array[11]
	;;	   ZLENS           DOUBLE    Array[11]
	;;	   SIGMA_PI_RATIO  DOUBLE    Array[11]
	;;	   RD              FLOAT     Array[11]
	;;	   RA              FLOAT     Array[11]
	;;	   H               FLOAT     Array[11]

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

  	chords={sigma_pi_ratio:sigma_pi_ratio,$	;;RATIO OF SIGMA LINES TO PI LINES  (0 IF NPA)
		 nchan:nchan,$				  		;;NUMBER OF CHANNELS
         chan_id:chan_id,$                  ;;CHANNEL ID (0 FOR FIDA, 1 FOR NPA)
		 xmid:xmid,$						;;X POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 ymid:ymid,$						;;Y POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
         zmid:zmid,$						;;Z POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 xlens:xlens,$						;;X POS. OF LENS/APERTURE [cm]
		 ylens:ylens,$						;;Y POS. OF LENS/APERTURE [cm]
 		 zlens:zlens,$						;;Z POS. OF LENS/APERTURE [cm]
  		 ra:ra,$      				        ;;RADIUS OF NPA APERTURE [cm] (0 IF FIDA)
  		 rd:rd,$      				        ;;RADIUS OF NPA DETECTOR [cm] (0 IF FIDA)
		 h:h}     		                    ;;SEPERATION BETWEEN APERTURE AND DETECTOR [cm] (0 IF FIDA)

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  vtor:vtor,$					;;TORODIAL ANGULAR VELOCITY [rad/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END 
```
As you can see it is all rather self-explanatory. All I need to do now to get ABCD to work with prefida is supply the above information.
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
