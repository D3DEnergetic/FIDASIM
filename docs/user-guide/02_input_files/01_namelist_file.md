title: Namelist File

#Namelist File
The namelist file contains all the basic settings. 
This file can be edited without having to re-create the other input files.

##Header
The header has some basic information about the simulation run. 
```fortran
!! Created: Thu Mar 24 04:59:13 2016
!! FIDASIM version: 1.0.0
!! Comment: This is a comment
```
Generally, "This is a comment" is not good comment.
A good comment should be descriptive.

##Shot Information
```fortran
&fidasim_inputs
shot =      1    !! Shot Number
time =  0.00051    !! Time [s]
runid = '051'    !! runID
result_dir = '/home/lstagner/'    !! Result Directory
```
These are pretty self-explainatory. I would note that the `runid` could be anything so .
In this particular example the `runid` was set to be `051` which is the time in ms.
This run was a part of a time series so choosing `051` made sorting the runs by time very easy.

##Input Files
```fortran
!! Input Files
tables_file = '/home/lstagner/FIDASIM/tables/atomic_tables.h5' !! Atomic Tables File
equilibrium_file = '/home/lstagner/051_equilibrium.h5'         !! File containing plasma parameters and fields
geometry_file = '/home/lstagner/051_geometry.h5'               !! File containing NBI and diagnostic geometry
distribution_file = '/home/lstagner/051_distribution.h5'       !! File containing fast-ion distribution
neutrals_file = '/home/lstagner/051_neutrals.h5'               !! File containing the neutral density
```
@note 
The neutrals file is actually an [output file](../03_output_files/index.html).
It is listed here since there is an option to load the neutral density from a file instead of calculating it. 

##Simulation Switches
```fortran
!! Simulation Switches
calc_bes =  1        !! Calculate Beam Emission and Halo Spectra
calc_brems =  1      !! Calculate Bremsstrahlung
calc_fida =  1       !! Calculate FIDA Spectra
calc_npa =  0        !! Calculate NPA
calc_birth =  1      !! Calculate Birth Profile
calc_fida_wght =  0  !! Calculate FIDA weights
calc_npa_wght =  0   !! Calculate NPA weights
load_neutrals =  0   !! Load neutrals from a preexisting neutrals file
dump_dcx =  1        !! Dump DCX neutrals and spectra
verbose =  1         !! Verbose
```

##Monte Carlo Settings
```fortran
!! Monte Carlo Settings
n_fida = 25000000 !! Number of FIDA mc particles
n_npa = 2500000   !! Number of NPA mc particles
n_nbi = 50000     !! Number of NBI mc particles
n_halo = 500000   !! Number of HALO mc particles
n_dcx = 500000    !! Number of DCX mc particles
n_birth = 10000   !! Number of BIRTH mc particles
```

##Neutral Beam Settings
```fortran
!! Neutral Beam Settings
ab =   1.00728             !! Beam Species mass [amu]
pinj =     1.326           !! Beam Power [MW]
einj =    15.000           !! Beam Energy [keV]
species_mix(1) =   0.85000 !! Beam Species Mix (Full component)
species_mix(2) =   0.10000 !! Beam Species Mix (Half component)
species_mix(3) =   0.05000 !! Beam Species Mix (Third component)
```

##Plasma Settings
```fortran
!! Plasma Settings
ai = 2.01411        !! Ion Species mass [amu]
impurity_charge = 6 !! Impurity Charge
```

##Beam Grid Settings
```fortran
!! Beam Grid Settings
nx =  80        !! Number of cells in X direction (Into Plasma)
ny =  80        !! Number of cells in Y direction
nz = 120        !! Number of cells in Z direction
xmin =  -80.000 !! Minimum X value [cm]
xmax =   80.000 !! Maximum X value [cm]
ymin =  -80.000 !! Minimum Y value [cm]
ymax =   80.000 !! Maximum Y value [cm]
zmin = -120.000 !! Minimum Z value [cm]
zmax =  120.000 !! Maximum Z value [cm]

!! Tait-Bryan Angles for z-y`-x`` rotation
alpha = 0.00000 !! Rotation about z-axis [rad]
beta  = 0.00000 !! Rotation about y`-axis [rad]
gamma = 0.00000 !! Rotation about x"-axis [rad]

!! Beam Grid origin in machine coordinates (cartesian)
origin(1) = 0.000 !! X value [cm]
origin(2) = 0.000 !! Y value [cm]
origin(3) = 0.000 !! Z value [cm]
```

##Wavelength Grid Settings
```fortran
!! Wavelength Grid Settings
nlambda = 2000      !! Number of Wavelengths
lambdamin = 647.000 !! Minimum Wavelength [nm]
lambdamax = 667.000 !! Maximum Wavelength [nm]
```

##Weight Function Settings
```fortran
!! Weight Function Settings
ne_wght = 50             !! Number of Energies for Weights
np_wght = 50             !! Number of Pitches for Weights
nphi_wght = 100          !! Number of Gyro-angles for Weights
emax_wght = 100.00       !! Maximum Energy for Weights [keV]
nlambda_wght = 1000      !! Number of Wavelengths for Weights 
lambdamin_wght = 647.000 !! Minimum Wavelength for Weights [nm]
lambdamax_wght = 667.000 !! Maximum Wavelength for Weights [nm]

/
```

#Example Namelist File
```fortran
!! Created: Thu Mar 24 04:59:13 2016
!! FIDASIM version: 1.0.0
!! Comment: This is a comment
&fidasim_inputs
shot =      1    !! Shot Number
time =  0.00051    !! Time [s]
runid = '051'    !! runID
result_dir = '/home/lstagner/'    !! Result Directory

!! Input Files
tables_file = '/home/lstagner/FIDASIM/tables/atomic_tables.h5' !! Atomic Tables File
equilibrium_file = '/home/lstagner/051_equilibrium.h5'         !! File containing plasma parameters and fields
geometry_file = '/home/lstagner/051_geometry.h5'               !! File containing NBI and diagnostic geometry
distribution_file = '/home/lstagner/051_distribution.h5'       !! File containing fast-ion distribution
neutrals_file = '/home/lstagner/051_neutrals.h5'               !! File containing the neutral density

!! Simulation Switches
calc_bes =  1        !! Calculate Beam Emission and Halo Spectra
calc_brems =  1      !! Calculate Bremsstrahlung
calc_fida =  1       !! Calculate FIDA Spectra
calc_npa =  0        !! Calculate NPA
calc_birth =  1      !! Calculate Birth Profile
calc_fida_wght =  0  !! Calculate FIDA weights
calc_npa_wght =  0   !! Calculate NPA weights
load_neutrals =  0   !! Load neutrals from a preexisting neutrals file
dump_dcx =  1        !! Dump DCX neutrals and spectra
verbose =  1         !! Verbose

!! Monte Carlo Settings
n_fida = 25000000 !! Number of FIDA mc particles
n_npa = 2500000   !! Number of NPA mc particles
n_nbi = 50000     !! Number of NBI mc particles
n_halo = 500000   !! Number of HALO mc particles
n_dcx = 500000    !! Number of DCX mc particles
n_birth = 10000   !! Number of BIRTH mc particles

!! Neutral Beam Settings
ab =   1.00728             !! Beam Species mass [amu]
pinj =     1.326           !! Beam Power [MW]
einj =    15.000           !! Beam Energy [keV]
species_mix(1) =   0.85000 !! Beam Species Mix (Full component)
species_mix(2) =   0.10000 !! Beam Species Mix (Half component)
species_mix(3) =   0.05000 !! Beam Species Mix (Third component)

!! Plasma Settings
ai = 2.01411        !! Ion Species mass [amu]
impurity_charge = 6 !! Impurity Charge

!! Beam Grid Settings
nx =  80        !! Number of cells in X direction (Into Plasma)
ny =  80        !! Number of cells in Y direction
nz = 120        !! Number of cells in Z direction
xmin =  -80.000 !! Minimum X value [cm]
xmax =   80.000 !! Maximum X value [cm]
ymin =  -80.000 !! Minimum Y value [cm]
ymax =   80.000 !! Maximum Y value [cm]
zmin = -120.000 !! Minimum Z value [cm]
zmax =  120.000 !! Maximum Z value [cm]

!! Tait-Bryan Angles for z-y`-x`` rotation
alpha = 0.00000 !! Rotation about z-axis [rad]
beta  = 0.00000 !! Rotation about y`-axis [rad]
gamma = 0.00000 !! Rotation about x"-axis [rad]

!! Beam Grid origin in machine coordinates (cartesian)
origin(1) = 0.000 !! X value [cm]
origin(2) = 0.000 !! Y value [cm]
origin(3) = 0.000 !! Z value [cm]

!! Wavelength Grid Settings
nlambda = 2000      !! Number of Wavelengths
lambdamin = 647.000 !! Minimum Wavelength [nm]
lambdamax = 667.000 !! Maximum Wavelength [nm]

!! Weight Function Settings
ne_wght = 50             !! Number of Energies for Weights
np_wght = 50             !! Number of Pitches for Weights
nphi_wght = 100          !! Number of Gyro-angles for Weights
emax_wght = 100.00       !! Maximum Energy for Weights [keV]
nlambda_wght = 1000      !! Number of Wavelengths for Weights 
lambdamin_wght = 647.000 !! Minimum Wavelength for Weights [nm]
lambdamax_wght = 667.000 !! Maximum Wavelength for Weights [nm]

/
```
