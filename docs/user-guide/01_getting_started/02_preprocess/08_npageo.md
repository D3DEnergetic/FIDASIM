title: NPA Geometry

#NPA Geometry Structure

[TOC]

#Structure Variables
**nchan**: Number of Channels

* type: `Int32`
* rank: 0

**system**: Name of NPA system(s)

* type: `String`
* rank: 0 or 1
* dims: [`nsystem`]

**data_source**: Source of NPA geometry

* type: `String`
* rank: 0

**a_shape**:Shape of the aperture

* type: `Int16`
* rank: 0

**d_shape**:Shape of the detector

* type: `Int16`
* rank: 0

The shapes of the detector and aperture can take the values

1. Rectangular
2. Circular

**a_cent**: Center of the aperture

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**a_redge**: Center of the apertures right edge

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**a_tedge**: Center of the apertures top edge

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**d_cent**: Center of the detector

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**d_redge**: Center of the detectors right edge

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**d_tedge**: Center of the detectors top edge

* type: `Float64`
* rank: 2
* dims: [3,`nchan`]
* units: cm

**radius**: Line of sight radius at midplane or tangency point

* type: `Float64`
* rank: 1
* dims: [`nchan`]
* units: cm

#Useful Routines

##check_npa.pro
[check_npa.pro](|url|/sourcefile/check_npa.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the NPA geometry structure has the correct format.

#Example Structure
```
IDL> help, npa
** Structure <48833a8>, 12 tags, length=512, data length=504, refs=1:
   NCHAN           LONG                 3
   SYSTEM          STRING    'NPA'
   DATA_SOURCE     STRING    'test_npa.pro'
   A_SHAPE         INT       Array[3]
   D_SHAPE         INT       Array[3]
   A_CENT          DOUBLE    Array[3, 3]
   A_REDGE         DOUBLE    Array[3, 3]
   A_TEDGE         DOUBLE    Array[3, 3]
   D_CENT          DOUBLE    Array[3, 3]
   D_REDGE         DOUBLE    Array[3, 3]
   D_TEDGE         DOUBLE    Array[3, 3]
   RADIUS          DOUBLE    Array[3]
```
