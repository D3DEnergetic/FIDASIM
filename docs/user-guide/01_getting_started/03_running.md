title: Running FIDASIM

#Running FIDASIM

[TOC]

Running FIDASIM is as easy as running

```
lstagner@computer:~/FIDASIM$ ./fidasim 
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ / 
/_/  /___//____//_/ |_|/___//___//_/  /_/  
                                           
Version: 1.0.0

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/

usage: ./fidasim namelist_file [num_threads]
```
Actually having FIDASIM produce something takes a bit of thought

##Recommended Hardware

The following settings will give a reasonable runtime.

* 32 threads on a shared memory node (All calculations are done on a single node)
* At least 2 GB of memory

@warning
By default FIDASIM will use all the threads available. If another process is hogging a core it will cause FIDASIM to stall. To prevent this use the `num_threads` optional argument

##Running Interactively

```
[lstagner@dawson061]% ./fidasim /p/fida/lstagner/TEST/test_1a_inputs.dat 16
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ / 
/_/  /___//____//_/ |_|/___//___//_/  /_/  
                                           
Version: 1.0.0

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/

---- Shot settings ----
 Shot:        1
 Time: 1000 [ms]
 Runid: test_1a
 
---- Input files ----
 Tables file: /u/lstagner/FIDASIM/tables/atomic_tables.h5
 Geometry file: /p/fida/lstagner/TEST/test_1a_geometry.h5
 Equilibrium file: /p/fida/lstagner/TEST/test_1a_equilibrium.h5
 Distribution file: /p/fida/lstagner/TEST/test_1a_distribution.h5
 
---- OpenMP settings ----
 Number of threads: 16
 
---- Beam grid settings ----
 Nx:  50
 Ny:  60
 Nz:  70
 dV:  8.00 [cm^3]
 alpha:  0.00 [rad]
 beta:   0.00 [rad]
 gamma:  0.00 [rad]
 origin: [   0.00,   0.00,   0.00] [cm]
 
---- Interpolation grid settings ----
 Nr:  70
 Nz: 100
 dA: 4.10 [cm^2]
 
---- Neutral beam settings ----
 Beam: test_beam               
 Power:    1.70 [MW]
 Voltage: 72.50 [keV]
 
---- Atomic tables settings ----
 Maximum n/m:  6
 Beam/Fast-ion mass:  2.014 [amu]
 Thermal/Bulk-ion mass:  2.014 [amu]
 Impurity mass: 12.011 [amu]
 
---- Fast-ion distribution settings ----
 Distribution type: Fast-ion Density Function F(energy,pitch,R,Z)
 Nenergy =   6
 Npitch  =   6
 Energy range = [67.33, 75.44]
 Pitch  range = [-0.10, 0.10]
 
---- FIDA/BES settings ----
 FIDA/BES System: SPECTRAL
 Number of channels:   3
 
---- NPA settings ----
 NPA System: NPA
 Number of channels:   3
 Calculating hit probabilities for NPA channels
                                                  
ndmc:     1:43:23
     # of markers:     50000
   birth profile written to: /p/fida/lstagner/TEST/test_1a_birth.h5
                              
dcx:      1:43:41
     # of markers:    505020
   dcx written to: /p/fida/lstagner/TEST/test_1a_dcx.h5
                              
halo:     1:44:32
     # of markers:    505180
     # of markers:    310573
     # of markers:    188148
     # of markers:    110872
     # of markers:     62806
     # of markers:     32484
     # of markers:     13881
   neutral density written to: /p/fida/lstagner/TEST/test_1a_neutrals.h5
                              
bremsstrahlung:     1:46:25
                              
fida:     1:46:25
     # of markers:   5049813
                              
   Spectra written to: /p/fida/lstagner/TEST/test_1a_spectra.h5

npa:     1:47:46
     # of markers:    505074
Number of NPA particles that hit a detector:   125638
                              
   NPA data written to: /p/fida/lstagner/TEST/test_1a_npa.h5

fida weight function:     1:49:46
 Number of Channels:   3
 Nlambda: 1000
 Nenergy:  50
 Npitch:  50
 Ngyro: 100
 Maximum Energy:  100.00
 LOS averaged: True
 
   Channel:   1
   Radius:  200.00
   Mean Fast-ion Density:    7.97429E+11
 
   Channel:   2
   Radius:  170.00
   Mean Fast-ion Density:    7.98346E+11
 
   Channel:   3
   Radius:  140.00
   Mean Fast-ion Density:    7.98330E+11
 
   FIDA weights written to: /p/fida/lstagner/TEST/test_1a_fida_weights.h5
                              
npa weight function:     1:50:02
 Number of Channels:   3
 Nenergy:  50
 Npitch:  50
 Maximum energy:  100.00
 
   Channel:   1
   Radius:    200.000
   Flux:      1.22243E+14
   Weight:    3.79893E+03
 
   Channel:   2
   Radius:    170.000
   Flux:      1.07364E+14
   Weight:    1.85565E+03
 
   Channel:   3
   Radius:    140.000
   Flux:      3.46488E+13
   Weight:    8.81099E+02
 
   NPA weights written to: /p/fida/lstagner/TEST/test_1a_npa_weights.h5
                              
END: hour, minute, second:  1:53:07
duration:                   0:15:53
```

##Submitting to a SGE scheduled cluster using `submit_fidasim`

`submit_fidasim` is a python routine that submits a FIDASIM job to a SGE cluster. For example

```
lstagner@computer:~$ submit_fidasim /u/lstagner/TEST

```
will submit any incomplete FIDASIM runs in the `/u/lstagner/TEST` directory. Alternatively
```

lstagner@computer:~$ submit_fidasim /u/lstagner/TEST/test_1a_inputs.dat

```
will submit just the `test_1a` FIDASIM run. See below for all the available options

```
lstagner@computer:~$ submit_fidasim -h
usage: submit_fidasim [-h] [-w WALLTIME] [-n NODES] [-ppn PPN] [-mem MEMORY]
                      [-ex EXECUTABLE] [-log LOG] [-pbs PBS] [-pre PRECALL]
                      [-post POSTCALL] [-rids RUNIDS [RUNIDS ...]] [-j JOBS]
                      [-c] [-v] [-db]
                      path

Creates a FIDASIM PBS job script and submits it using qsub

positional arguments:
  path                  Namelist file or result directory

optional arguments:
  -h, --help            show this help message and exit
  -w WALLTIME, --walltime WALLTIME
                        Set walltime. Defaults to 5:00:00
  -n NODES, --nodes NODES
                        Set number of nodes. Defaults to 1
  -ppn PPN              Set processors per node. Defaults to 8
  -mem MEMORY, --memory MEMORY
                        Set required memory. Defaults to 2048mb
  -ex EXECUTABLE, --executable EXECUTABLE
                        Set path to FIDASIM executable. Defaults to
                        /home/lstagner/FIDASIM/fidasim
  -log LOG              Set log directory. Defaults to result directory.
  -pbs PBS              Additional PBS directive
  -pre PRECALL, --precall PRECALL
                        Command to run before code execution
  -post POSTCALL, --postcall POSTCALL
                        Command to run after code execution
  -rids RUNIDS [RUNIDS ...], --runids RUNIDS [RUNIDS ...]
                        List of run ids, accepts regex
  -j JOBS, --jobs JOBS  Split runs into N jobs. Defaults to 1 job per namelist
  -c, --clobber         Overwrite existing runs
  -v, --verbose         Verbose
  -db, --debug          Debug mode. Does not submit job
```
