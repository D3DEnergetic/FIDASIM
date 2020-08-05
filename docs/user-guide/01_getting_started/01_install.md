title: Installation

#Installing FIDASIM
So you have decided to install FIDASIM. Don't worry this should be relatively painless.

@note
If you experiance problems installing FIDASIM you did something wrong and you should [**let us know**](https://github.com/D3DEnergetic/FIDASIM/issues/new) so we can laugh at you (and also help you)

The following code snippets assumes you are using a BASH shell.
To figure out what shell you currently have run `echo $SHELL` to find out.

[TOC]

---

##Dependencies
FIDASIM has the following dependencies:

* Linux because all other operating systems are inferior
* [Make](https://www.gnu.org/software/make/) for building FIDASIM. (Required)
* Fortran, C, and C++ compilers. (Required)
  [GNU(>v5)](https://gcc.gnu.org/) and [Intel(>13)](https://software.intel.com/en-us/intel-compilers) compilers are supported.
  Note you cannot mix and match different compilers.
* [zlib](http://zlib.net/) for file compression. (Required)
* [Anaconda Python](https://www.continuum.io/why-anaconda) for python scripts and pre-processing. (Optional)
* [IDL](http://www.harrisgeospatial.com/ProductsandSolutions/GeospatialProducts/IDL.aspx) for pre-processing (Optional)
* [HDF5 1.8.16](https://www.hdfgroup.org/HDF5/) for reading and writing data files (Included)
* [EFIT](https://fusion.gat.com/theory/Efit) for electro-magnetic fields (Partially Included)
* [git](https://git-scm.com/) for version control (Optional)
* [Ford](https://github.com/cmacmackin/ford) for creating HTML documentation (Optional)
* [LinkChecker](http://wummel.github.io/linkchecker/) for validating HTML documentation (Optional)

The following commands will install the required dependencies (Tested on Ubuntu 16.04)
```bash
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install gfortran
sudo apt-get install zlib1g-dev
```

##Getting FIDASIM source
It's rather difficult to run software you haven't downloaded. There are two ways of getting the source code.

###Downloading source directly
The most recent version of FIDASIM ({!../VERSION!}) can be downloaded from [here](https://github.com/D3DEnergetic/FIDASIM/releases)
Once you have downloaded the `.tar.gz` or `.zip` file unpack it using the following commands.
```bash
tar -zxf FIDASIM-{!../VERSION!}.tar.gz
```
or if you downloaded the `.zip` file
```bash
unzip FIDASIM-{!../VERSION!}.zip
```

There should now be a directory named `FIDASIM-{!../VERSION!}` in your current directory. Enter the directory using
```bash
cd FIDASIM-{!../VERSION!}
```

###Cloning the git repository
If you are planning to develop FIDASIM it is best to use git.
```bash
git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
cd FIDASIM
```

By default you will be on the master branch which may or may not be in a usable state.
To use the release version run the command
```bash
git checkout v{!../VERSION!}
```

##Setting up environmental variables
FIDASIM needs to know where some things are so you have to set the following environmental variables in your `.bashrc` file located in your home directory.
```bash
export FIDASIM_DIR=/path/to/fidasim/install/directory
export FC=gfortran #use 'ifort' for Intel Fortran compiler
export CC=gcc      #use 'icc' for Intel C compiler
export CXX=g++     #use 'icpc' for Intel C++ compiler

#For using helper routines
export PATH=$FIDASIM_DIR/deps/hdf5/bin:$FIDASIM_DIR/lib/scripts:$PATH
export PATH=$FIDASIM_DIR/deps/efit:$FIDASIM_DIR/test:$PATH
export IDL_PATH="+$FIDASIM_DIR:$IDL_PATH:<IDL_DEFAULT>"
export PYTHONPATH=$FIDASIM_DIR/lib/python:$PYTHONPATH

ulimit -s unlimited #Omit this if you like segfaults
```
replacing `/path/to/fidasim/install/directory` with the real directory. To set the environmental variables in the current shell run
```bash
source ~/.bashrc
```

##Building FIDASIM
Once you are in the source directory (and have all the dependencies installed) run the following
```bash
make
```
This will build the fidasim executable using the default compiler options. Run `make help` to view all the available options.
Once running, go get a coffee since it will take a while because HDF5 is being built as well.

Once make has completed check if FIDASIM compiled correctly.
```
user@computer:~/FIDASIM-{!../VERSION!}$ ./fidasim
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ /
/_/  /___//____//_/ |_|/___//___//_/  /_/

Version: {!../VERSION!}

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/
```
Good job! You installed FIDASIM! But wait theres more.

##Generating Atomic Tables
Calculating reaction rates on the fly is time consuming so we need to pre-compute them to save time.
The following code snippit will generate the atomic tables using the default settings.
The default settings should be appropriate for most use cases, however, it may be necessary to generate custom atomic tables.
In that case edit the file `tables/table_settings.dat` before running the following command if FIDASIM was compiled with OpenMP (the default build)
```bash
./tables/generate_tables ./tables/default_settings.dat [num_threads]
```
or if FIDASIM was built with MPI
```bash
mpirun -np num_processes ./tables/generate_tables ./tables/default_settings.dat
```
@warning
This is computationally expensive so make sure you run this on a computer
where you won't get angry emails for using up all the CPU's

Now would be a good time to get more coffee... or maybe a nap.

##Run a test case
From the command line
```bash
run_tests.py, "/place/where/you/want/the/output"
```
** Note: This requires python **

Or from within IDL
```idl
IDL> run_tests, "/place/where/you/want/the/output"
```

Some stuff that will make sense later will flash by and when its done you should see something like
```text
SUCCESS: FIDASIM pre-processing completed
To run FIDASIM use the following command
/u/lstagner/FIDASIM/fidasim /p/fida/lstagner/TEST/test_1a_inputs.dat

```

Now do what the computer says.
Think of as good practice for when the [robots take over](https://www.youtube.com/watch?v=7Pq-S557XQU).

It should print out the following.
```
[lstagner@dawson061]% /u/lstagner/FIDASIM/fidasim /p/fida/lstagner/TEST/test_1a_inputs.dat
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ /
/_/  /___//____//_/ |_|/___//___//_/  /_/

Version: {!../VERSION!}

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/

---- Shot settings ----
 Shot:        1
 Time: 1000 [ms]
 Runid: test

---- Input files ----
 Tables file: /p/fida/FIDASIM/tables/atomic_tables.h5
 Geometry file: /p/fida/lstagner/TEST/test_geometry.h5
 Equilibrium file: /p/fida/lstagner/TEST/test_equilibrium.h5
 Distribution file: /p/fida/lstagner/TEST/test_distribution.h5

---- OpenMP settings ----
 Number of threads: 16

---- Atomic tables settings ----
 Maximum n/m:  6
 Beam/Fast-ion mass:  2.014 [amu]
 Thermal/Bulk-ion mass:  2.014 [amu]
 Impurity mass: 12.011 [amu]

---- Interpolation grid settings ----
 Nr:  70
 Nz: 100
 Nphi:   1
 dA:  4.00 [cm^2]

---- Beam grid settings ----
 Nx:  50
 Ny:  60
 Nz:  70
 dV:  8.00 [cm^3]
 alpha:  0.00 [rad]
 beta:   0.00 [rad]
 gamma:  0.00 [rad]
 origin: [   0.00,   0.00,   0.00] [cm]
 Number of cells in plasma:   184494

---- Neutral beam settings ----
 Beam: test_beam
 Power:    1.70 [MW]
 Voltage:  72.50 [keV]

---- Passive grid settings ----
 Nr:  70
 Nz: 100
 Nphi:  10
 R  range = [100.00,238.00]
 Z  range = [-100.00, 98.00]
 Phi  range = [ 4.25, 5.15]
 dA:  4.00 [cm^3]

---- Fast-ion distribution settings ----
 Distribution type: Fast-ion Density Function F(energy,pitch,R,Z,Phi)
 Nenergy =  75
 Npitch  =  50
 Nr  =  70
 Nz  = 100
 Nphi  =   1
 Energy range = [ 0.81,120.87]
 Pitch  range = [-0.98, 0.98]
 R  range = [100.00,238.00]
 Z  range = [-100.00, 98.00]
 Phi  range = [ 0.00, 0.00]
 Ntotal =  1.214E+19

---- FIDA/BES settings ----
 FIDA/BES System: SPECTRAL
 Number of channels:     3

---- NPA settings ----
 NPA System: NPA
 Number of channels:   3

 nbi:     19:06:28 --- elapsed: 0:00:23
     # of markers:     50000
   birth profile written to: /p/fida/lstagner/TEST/test_birth.h5

 dcx:     19:06:42 --- elapsed: 0:00:37
     # of markers:    500000

 halo:    19:07:13 --- elapsed: 0:01:08
     # of markers:    527175 --- Seed/DCX: 1.000
     # of markers:    672040 --- Seed/DCX: 0.600
     # of markers:    792640 --- Seed/DCX: 0.366
     # of markers:    873160 --- Seed/DCX: 0.225
     # of markers:    900615 --- Seed/DCX: 0.139
     # of markers:    913130 --- Seed/DCX: 0.086
     # of markers:    919085 --- Seed/DCX: 0.053
     # of markers:    921060 --- Seed/DCX: 0.033
     # of markers:    921935 --- Seed/DCX: 0.020
     # of markers:    922180 --- Seed/DCX: 0.013

 write neutrals:    19:14:31 --- elapsed: 0:08:26
   neutral density written to: /p/fida/lstagner/TEST/test_neutrals.h5

 cold:    19:14:35 --- elapsed: 0:08:30

 bremsstrahlung:    19:14:36 --- elapsed: 0:08:31

 fida:    19:14:36 --- elapsed: 0:08:31
     # of markers:   5000000

 pfida:   19:15:19 --- elapsed: 0:09:14
     # of markers:  50000000

 write spectra:    19:20:47 --- elapsed: 0:14:42
   Spectra written to: /p/fida/lstagner/TEST/test_spectra.h5

 npa:     19:20:47 --- elapsed: 0:14:42
     # of markers:      5000000
   Number of Active NPA particles that hit a detector:    40733

 pnpa:     19:21:21 --- elapsed: 0:15:16
     # of markers:     50000000
   Number of Passive NPA particles that hit a detector:   116510

 write npa:    19:26:47 --- elapsed: 0:20:42
   NPA data written to: /p/fida/lstagner/TEST/test_npa.h5

 neutron rate:    19:26:47 --- elapsed: 0:20:42
   Rate:      5.97592E+14 [neutrons/s]

 write neutrons:    19:28:16 --- elapsed: 0:22:11
   Neutrons written to: /p/fida/lstagner/TEST/test_neutrons.h5

 fida weight function:    19:28:25 --- elapsed: 0:22:20
  Number of Channels:     3
  Nlambda: 1000
  Nenergy:  50
  Npitch:  50
  Ngyro: 100
  Maximum Energy:  100.00
  LOS averaged: True

   Channel:     1
   Radius:  200.00
   Mean Fast-ion Density:    5.00757E+11

   Channel:     2
   Radius:  170.00
   Mean Fast-ion Density:    5.00757E+11

   Channel:     3
   Radius:  140.00
   Mean Fast-ion Density:    5.00757E+11

 write fida weights:    19:28:41 --- elapsed: 0:22:36
   FIDA weights written to: /p/fida/lstagner/TEST/test_fida_weights.h5

 npa weight function:    19:28:48 --- elapsed: 0:22:43
  Number of Channels:   3
  Nenergy:  50
  Npitch:  50
  Maximum energy:  100.00

   Channel:   1
   Radius:    200.000
   Flux:      0.00000E+00
   Weight:    3.48351E+03

   Channel:   2
   Radius:    170.000
   Flux:      0.00000E+00
   Weight:    1.41572E+03

   Channel:   3
   Radius:    140.000
   Flux:      0.00000E+00
   Weight:    8.38209E+02

 write npa weights:    19:33:50 --- elapsed: 0:27:45
   NPA weights written to: /p/fida/lstagner/TEST/test_npa_weights.h5

 END: hour:minute:second 19:33:50 --- elapsed: 0:27:45
```

Congratulations! You followed the instructions.

##Now what
Most likely you wont be satisfied by just running a test case. Click [here](./02_preprocess.html) to learn how to make the input files used by FIDASIM.
