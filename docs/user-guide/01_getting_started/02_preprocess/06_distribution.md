title: Fast-ion Distribution

#Distribution Structure

[TOC]

#Useful Routines

##extract_transp_fbm
`extract_transp_fbm` is a python routine that extracts a fast-ion distribution function from a TRANSP/NUBEAM run.
```
lstagner@computer:~/FIDASIM$ extract_transp_fbm -h
usage: extract_transp_fbm [-h] [-fid FILE_ID [FILE_ID ...]] [-d DEVICE]
                          [-o OUTPUT_DIR] [-t TYPE] [-v] [-db]
                          path runid

Extracts fast-ion distribution CDF file from a TRANSP run

positional arguments:
  path                  Result directory
  runid                 TRANSP run ID

optional arguments:
  -h, --help            show this help message and exit
  -fid FILE_ID [FILE_ID ...], --file_id FILE_ID [FILE_ID ...]
                        File ID list e.g. 1 2 3 ...
  -d DEVICE, --device DEVICE
                        Set device
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default = cwd
  -t TYPE, --type TYPE  Type of distribution: c: Guiding Center (default) or
                        p: Particle Location
  -v, --verbose         Verbose
  -db, --debug          Debug mode
```

##read_nubeam.pro
[read_nubeam.pro](|url|/sourcefile/read_nubeam.pro.html) is an IDL routine that reads the TRANSP/NUBEAM distribution file produced by the `get_fbm` routine.

##read_mc_nubeam.pro
[read_mc_nubeam.pro](|url|/sourcefile/read_mc_nubeam.pro.html) is an IDL routine reads the TRANSP/NUBEAM Guiding Center Monte Carlo distribution file produced by the `get_fbm` routine.

##read_spiral.pro
[read_spiral.pro](|url|/sourcefile/read_spiral.pro.html) is an IDL routine reads the Guiding Center Monte Carlo distribution file produced by SPIRAL.

##check_distribution.pro
[check_distribution.pro](|url|/sourcefile/check_distribution.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the distribution structure has the correct format.

#Example Structure

##Distribution Function
```
IDL> help, dist
** Structure <44833f8>, 9 tags, length=2072144, data length=2072126, refs=1:
   TYPE            INT              1
   TIME            DOUBLE           1.0000000
   NENERGY         INT              6
   ENERGY          DOUBLE    Array[6]
   NPITCH          INT              6
   PITCH           DOUBLE    Array[6]
   F               DOUBLE    Array[6, 6, 70, 100]
   DENF            DOUBLE    Array[70, 100]
   DATA_SOURCE     STRING    '/u/lstagner/FIDASIM/test/test_fi_1.cdf'
```
##Guiding Center Monte Carlo Distribution
```
IDL> help,a
** Structure <1da63648>, 11 tags, length=210000040, data length=210000032, refs=1:
   TYPE            INT              2
   TIME            DOUBLE          0.81199998
   DATA_SOURCE     STRING    '/u/lstagner/mc_159243H06_9'
   NPARTICLE       LONG           5000000
   NCLASS          INT              1
   R               DOUBLE    Array[5000000]
   Z               DOUBLE    Array[5000000]
   ENERGY          DOUBLE    Array[5000000]
   PITCH           DOUBLE    Array[5000000]
   CLASS           INT       Array[5000000]
   WEIGHT          DOUBLE    Array[5000000]
```

##Full-Orbit Monte Carlo Distribution
```
IDL> help,a
** Structure <f10daa8>, 12 tags, length=9516944, data length=9516932, refs=1:
   TYPE            INT              3
   DATA_SOURCE     STRING    'full_orbit.h5'
   TIME            DOUBLE        0.0030000024
   NPARTICLE       LONG            190338
   NCLASS          INT              1
   R               DOUBLE    Array[190338]
   Z               DOUBLE    Array[190338]
   VR              DOUBLE    Array[190338]
   VT              DOUBLE    Array[190338]
   VZ              DOUBLE    Array[190338]
   CLASS           INT       Array[190338]
   WEIGHT          DOUBLE    Array[190338]
```

#Structure Variables
**time**: Distribution time

* type: `Float64`
* rank: 0
* units: s

**data_source**: Source of fast-ion distribution data

* type: `String`
* rank: 0

**type**: Distribution Type

* type: `Int16`
* rank: 0

The variable `type` can take 3 different values

1. Fast-ion Distribution Function
2. Guiding Center Monte Carlo Distribution
3. Full-Orbit Monte Carlo Distribution

Depeding on the distribution `type` the remain structure variables can differ.

##Fast-ion Distribution Function
**nenergy**: Number of energy values

* type: `Int16`
* rank: 0

**npitch**: Number of pitch values

* type: `Int16`
* rank: 0

**energy**: Array of energy values

* type: `Float64`
* rank: 1
* dims: [`nenergy`]
* units: keV

**pitch**: Array of pitch values w.r.t. the magnetic field

* type: `Float64`
* rank: 1
* dims: [`npitch`]

**denf**: Fast-ion density

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm^-3

**f**: Fast-ion distribution function: F(E,p,R,Z)

* type: `Float64`
* rank: 4
* dims: [`nenergy`,`npitch`,`nr`,`nz`]
* units: fast-ions/(dE dP cm^3)

##Monte Carlo Distribution
**nparticle**: Number of MC particles

* type: `Int16`
* rank: 0

**nclass**: Number of orbit classes

* type: `Int16`
* rank: 0

**r**: R positions of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: cm

**z**: Z positions of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: cm

**weight**: Weights of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]

The sum(weight) = # of Fast-ions in phase space sampled by the MC particles

**class**: Orbit classes of the MC particles

* type: `Int16`
* rank: 1
* dims: [`nparticle`]

The class can take values in the range of 1:`nclass`

###Guiding Center
**energy**: Energy of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: keV

**pitch**: Pitch of the MC particles w.r.t. the magnetic field

* type: `Float64`
* rank: 1
* dims: [`nparticle`]

###Full-Orbit
**vr**: Radial velocity of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: cm/s

**vt**: Torodial velocity of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: cm/s

**vz**: Z velocity of the MC particles

* type: `Float64`
* rank: 1
* dims: [`nparticle`]
* units: cm/s

