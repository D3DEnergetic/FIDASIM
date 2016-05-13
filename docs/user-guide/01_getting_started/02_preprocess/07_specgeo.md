title: Spectral Geometry

#Spectral Geometry Structure

[TOC]

#Useful Routines

##check_spec.pro
[check_spec.pro](|url|/sourcefile/check_spec.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the spectral geometry structure has the correct format.

#Example Structure
```
IDL> help, spec
** Structure <4571cb8>, 8 tags, length=256, data length=252, refs=1:
   NCHAN           LONG                 3
   SYSTEM          STRING    'SPECTRAL'
   DATA_SOURCE     STRING    'test_chords.pro'
   ID              STRING    Array[3]
   LENS            DOUBLE    Array[3, 3]
   AXIS            DOUBLE    Array[3, 3]
   SPOT_SIZE       DOUBLE    Array[3]
   SIGMA_PI        DOUBLE    Array[3]
   RADIUS          DOUBLE    Array[3]
```

#Structure Variables
**nchan**: Number of channels

* type: `Int32`
* rank: 0

**system**: Name of spectroscopic system(s)

* type: `String`
* rank: 0 

**data_source**: Source of spectral geometry

* type: `String`
* rank: 0

**id**: Line of sight ID

* type: `String`
* rank: 1
* dims: [`nchan`]

**lens**: Lens location in machine coordinates

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**axis**: Optical axis of the lines of sight

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]

**spot_size**: Radius of the collecting volume

* type: `Float64`
* rank: 1
* dims: [`nchan`]
* units: cm

**sigma_pi**: Ratio of the intensities of the sigma and pi stark lines

* type: `Float64`
* rank: 1
* dims: [`nchan`]

**radius**: Line of sight radius at midplane or tangency point

* type: `Float64`
* rank: 1
* dims: [`nchan`]
* units: cm

