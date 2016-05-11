title: General Settings

#Inputs Structure
The `inputs` structure contains the basic information needed by FIDASIM

[TOC]

#Structure Variables

The variables can be split up into several sections

##Shot Information
**shot**: Shot Number

* type: `Int32`
* rank:  0

**time**: Shot time

* type: `Float64`
* rank: 0
* units: s

**runid**: FIDASIM run ID

* type: `String`
* rank: 0

**comment**: FIDASIM run comment 

* type: `String`
* rank: 0

##Directories and Files
**install_dir**: FIDASIM install directory

* type: `String`
* rank: 0

**result_dir**: Result directory

* type: `String`
* rank: 0

**tables_file**: Atomic Tables file

* type: `String`
* rank: 0

##Simulation Switches
**calc_bes**: Calculate Beam Emission and Halo spectra

* type: `Int16`
* rank: 0

**calc_brems**: Calculate Bremsstrahlung

* type: `Int16`
* rank: 0

**calc_fida**: Calculate FIDA Spectra

* type: `Int16`
* rank: 0

**calc_npa**: Calculate NPA

* type: `Int16`
* rank: 0

**calc_birth**: Calculate Birth Profile

* type: `Int16`
* rank: 0

**calc_fida_wght**: Calculate FIDA weights

* type: `Int16`
* rank: 0

**calc_npa_wght**: Calculate NPA weights

* type: `Int16`
* rank: 0

**load_neutrals**: Load neutrals from a preexisting neutrals file

* type: `Int16`
* rank: 0

**dump_dcx**: Dump DCX neutrals and spectra

* type: `Int16`
* rank: 0

**verbose**: Verbose output

* type: `Int16`
* rank: 0

The simulation switches can take values 0, 1, or 2. A value of zero and one will turn the calculation off and on respectively. 
A value of two will turn on additional functionality.

##Monte Carlo Settings
**n_fida**: Number of FIDA mc particles

* type: `Int32`
* rank: 0

**n_npa**: Number of NPA mc particles

* type: `Int32`
* rank: 0

**n_nbi**: Number of NBI mc particles 

* type: `Int32`
* rank: 0

**n_halo**: Number of HALO mc particles 

* type: `Int32`
* rank: 0

**n_dcx**: Number of DCX mc particles

* type: `Int32`
* rank: 0

**n_birth**: Number of BIRTH mc particles

* type: `Int32`
* rank: 0

##Neutral Beam Settings

**ab**: Beam Species mass

* type: `Float64`
* rank: 0
* units: amu

**pinj**: Beam Power

* type: `Float64`
* rank: 0
* units: MW

**einj**: Beam Energy

* type: `Float64`
* rank: 0
* units: keV

**species_mix**: Beam Species Mix / Current Fractions (Full, Half, Third components)

* type: `Float64`
* rank: 1
* dims: [3]

##Plasma Settings
**ai**: Thermal Ion Species mass

* type: `Float64`
* rank: 0
* units: amu

**impurity_charge**: Impurity Charge

* type: `Int16`
* rank: 0

##Beam Grid Settings
**nx**: Number of cells in the X direction

* type: `Int16`
* rank: 0

**ny**: Number of cells in the Y direction

* type: `Int16`
* rank: 0

**nz**: Number of cells in the Z direction

* type: `Int16`
* rank: 0

**xmin**: Minimum X value

* type: `Float64`
* rank: 0
* units: cm

**xmax**: Maximum X value

* type: `Float64`
* rank: 0
* units: cm

**ymin**: Minimum Y value

* type: `Float64`
* rank: 0
* units: cm

**ymax**: Maximum Y value

* type: `Float64`
* rank: 0
* units: cm

**zmin**: Minimum Z value

* type: `Float64`
* rank: 0
* units: cm

**zmax**: Maximum Z value

* type: `Float64`
* rank: 0
* units: cm

**alpha**: Rotation angle about z-axis

* type: `Float64`
* rank: 0
* units: radians

**beta**: Rotation angle about the y'-axis

* type: `Float64`
* rank: 0
* units: radians

**gamma**: Rotation angle about the x"-axis

* type: `Float64`
* rank: 0
* units: radians

**origin**: Origin in machine coordinates

* type: `Float64`
* rank: 1
* dims: [3]
* units: cm

The coordinates system used depends on the values of the rotation angles (`alpha`, `beta`, `gamma`) and `origin` variables.
If all the rotation angles and `origin` are set to zero then the coordinate system is identical to machine coordinates. 
The angles variables and `origin` define a rotated coordinate system. 
The rotation angles and origin can best be described by an example. 

1. With your right hand point your index finger pointing in the +x direction with your middle finger and thumb pointing in the +y and +z direction respectively.
2. Rotate about your thumb (z-axis) by `alpha` (ccw = +angle, cw = -angle)
3. Rotate about your middle finger (y'-axis) by `beta`
4. Rotate about your index finger (x"-axis) by `gamma`
5. Move your right hand to the `origin`
6. Define `(x|y|z)_(min|max)` by this coordinate system with your index finger being the new +x-axis

##Wavelength Grid Settings
**nlambda**: Number of wavelengths

* type: `Int16`
* rank: 0

**lambdamin**: Minimum Wavelength

* type: `Float64`
* rank: 0
* units: nm

**lambdamax**: Maximum Wavelength

* type: `Float64`
* rank: 0
* units: nm

##Weight Function Settings
**ne_wght**: Number of energies

* type: `Int16`
* rank: 0

**np_wght**: Number of pitches

* type: `Int16`
* rank: 0

**nphi_wght**: Number of gyro-angles

* type: `Int16`
* rank: 0

**emax_wght**: Maximum energy

* type: `Float64`
* rank: 0
* units: keV

**nlambda_wght**: Number of wavelengths

* type: `Int16`
* rank: 0

**lambdamin_wght**: Minimum wavelength

* type: `Float64`
* rank: 0
* units: nm

**lambdamax_wght**: Maximum wavelength

* type: `Float64`
* rank: 0
* units: nm

#Useful Routines
##beam_grid.pro
It is convienient to define the grid to be aligned with the beam sightline. To faciliate this an IDL routine [beam_grid.pro](|url|/sourcefile/beam_grid.pro.html) is available to automatically calculate the beam-aligned grid definition.

##vars_to_struct.pro
The IDL routine [vars_to_struct.pro](|url|/sourcefile/vars_to_struct.pro.html) takes variables defined in the local scope and stores them into a structure. This is useful for creating input files.

##check_inputs.pro
[check_inputs.pro](|url|/sourcefile/check_inputs.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the inputs structure has the correct format.

#Example Structure
```
IDL> help, inputs
** Structure <44a4468>, 53 tags, length=392, data length=362, refs=3:
   SHOT            LONG                 1
   TIME            DOUBLE           1.0000000
   RUNID           STRING    'test_1a'
   DEVICE          STRING    'TEST'
   COMMENT         STRING    'Non-rotated, Non-tilted grid; flat profiles'
   INSTALL_DIR     STRING    '/home/lstagner/FIDASIM'
   RESULT_DIR      STRING    '/home/lstagner/TEST'
   TABLES_FILE     STRING    '/home/lstagner/FIDASIM/tables/atomic_tables.h5'
   NX              INT             50
   NY              INT             60
   NZ              INT             70
   XMIN            DOUBLE          -50.000000
   XMAX            DOUBLE           50.000000
   YMIN            DOUBLE          -230.00000
   YMAX            DOUBLE          -110.00000
   ZMIN            DOUBLE          -70.000000
   ZMAX            DOUBLE           70.000000
   ALPHA           DOUBLE           0.0000000
   BETA            DOUBLE           0.0000000
   GAMMA           DOUBLE           0.0000000
   ORIGIN          DOUBLE    Array[3]
   EINJ            DOUBLE           72.500000
   PINJ            DOUBLE           1.7000000               
   SPECIES_MIX     DOUBLE    Array[3]
   AB              DOUBLE           2.0141018
   AI              DOUBLE           2.0141078
   IMPURITY_CHARGE INT              6
   LAMBDAMIN       DOUBLE           647.00000
   LAMBDAMAX       DOUBLE           667.00000
   NLAMBDA         INT           2000
   N_FIDA          LONG           5000000
   N_NPA           LONG            500000
   N_NBI           LONG             50000
   N_HALO          LONG            500000
   N_DCX           LONG            500000
   N_BIRTH         LONG             10000
   NE_WGHT         INT             50
   NP_WGHT         INT             50
   NPHI_WGHT       INT            100
   EMAX_WGHT       DOUBLE           100.00000
   NLAMBDA_WGHT    INT           1000
   LAMBDAMIN_WGHT  DOUBLE           647.00000
   LAMBDAMAX_WGHT  DOUBLE           667.00000               
   CALC_NPA        INT              1
   CALC_BREMS      INT              1
   CALC_BES        INT              1
   CALC_FIDA       INT              1
   CALC_BIRTH      INT              1
   CALC_FIDA_WGHT  INT              1
   CALC_NPA_WGHT   INT              1
   LOAD_NEUTRALS   INT              0
   DUMP_DCX        INT              1
   VERBOSE         INT              1
```

