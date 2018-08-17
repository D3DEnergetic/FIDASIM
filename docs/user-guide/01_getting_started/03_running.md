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
                                           
Version: v2.0.0-dev-96-g30a264b

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/
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
[lstagner@dawson061]% ./fidasim /p/fida/lstagner/TEST/test_inputs.dat 16
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ / 
/_/  /___//____//_/ |_|/___//___//_/  /_/  
                                           
Version: v2.0.0-dev-96-g30a264b

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

##Submitting to a clusters job schedular using `submit_fidasim`

`submit_fidasim` is a python routine that schedules a FIDASIM job on a cluster. For example

```
lstagner@computer:~$ submit_fidasim /u/lstagner/TEST

```
will submit any incomplete FIDASIM runs in the `/u/lstagner/TEST` directory. Alternatively
```

lstagner@computer:~$ submit_fidasim /u/lstagner/TEST/test_inputs.dat

```
will submit just the `test` FIDASIM run.
Slurm and PBS resource managers are supported. Run `submit_fidasim -h` for the full documentation.
