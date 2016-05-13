title: Neutral Beam Geometry

#Neutral Beam Geometry Structure

[TOC]

#Useful Routines
## nubeam_geometry.pro
For convienience, the IDL routine [nubeam_geometry.pro](|url|/sourcefile/nubeam_geometry.pro.html) is available to transform the nubeam beam geometry definition.

##check_beam.pro
[check_beam.pro](|url|/sourcefile/check_beam.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the beam structure has the correct format.

#Example Structure
```
IDL> help, nbi
** Structure <4882ea8>, 11 tags, length=168, data length=162, refs=1:
   NAME            STRING    'test_beam'
   SHAPE           INT              1
   DATA_SOURCE     STRING    'test_beam.pro'
   SRC             DOUBLE    Array[3]
   AXIS            DOUBLE    Array[3]
   WIDY            DOUBLE           6.0000000
   WIDZ            DOUBLE           24.000000
   DIVY            DOUBLE    Array[3]
   DIVZ            DOUBLE    Array[3]
   FOCY            DOUBLE           999999.90
   FOCZ            DOUBLE           1000.0000
```
#Structure Variables
**name**: Beam Name

* type: `String`
* rank: 0

**shape**: Source Grid Shape

* type: `Int16`
* rank: 0

The shape of the source grid take the value of 

1. Rectangular
2. Circular

**data_source**: Geometry source

* type: `String`
* rank: 0

**src**: Position of source

* type: `Float64`
* rank: 1
* dims: [3]
* units: cm

**axis**: Beam sightline axis/direction

* type: `Float64`
* rank: 1
* dims: [3]

**widy**: Source grid half-width in the horizontal direction

* type: `Float64`
* rank: 0
* units: cm

**widz**: Source grid half-height in the vertical direction

* type: `Float64`
* rank: 0
* units: cm

**divy**: Horizontal beam divergence

* type: `Float64`
* rank: 1
* dims: [3]
* units: radians

**divz**: Vertical beam divergence

* type: `Float64`
* rank: 1
* dims: [3]
* units: radians

**focy**: Horizontal focal length

* type: `Float64`
* rank: 0
* units: cm

**focz**: Vertical focal length

* type: `Float64`
* rank: 0
* units: cm

