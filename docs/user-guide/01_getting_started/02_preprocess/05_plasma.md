title: Plasma Parameters

#Plasma Structure

[TOC]

#Structure Variables
**time**: Time when plasma parameters were extracted

* type: `Float64`
* rank: 0
* units: Seconds [s]

**data_source**: Plasma profiles data source

* type: `String`
* rank: 0

**mask**: Boolean mask that indicates where plasma parameters are known

* type: `Int16`
* rank: 2
* dims: [`nr`,`nz`]

**te**: Electron temperature

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: Kiloelecton-volt [keV]

**ti**: Ion temperature

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: Kiloelecton-volt [keV]

**vr**: Bulk plasma flow in the radial-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm/s

**vt**: Bulk plasma flow in the torodial/phi-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm/s

**vz**: Bulk plasma flow in the z-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm/s

**dene**: Electron number density

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm^-3

**zeff**: Effective nuclear charge

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]

#Useful Routines

##extract_transp_plasma.pro
The IDL routine [extract_transp_plasma.pro](|url|/sourcefile/extract_transp_plasma.pro.html) is provided to extract the plasma parameters structure from a TRANSP run. 

##check_plasma.pro
[check_plasma.pro](|url|/sourcefile/check_plasma.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the plasma structure has the correct format.

#Example Structure
```
IDL> help, plasma
** Structure <448c508>, 10 tags, length=406024, data length=406024, refs=1:
   TIME            DOUBLE           1.0000000
   DATA_SOURCE     STRING    '/home/lstagner/FIDASIM/test/test_profiles.pro'
   MASK            INT       Array[70, 100]
   TE              DOUBLE    Array[70, 100]
   TI              DOUBLE    Array[70, 100]
   VR              DOUBLE    Array[70, 100]
   VT              DOUBLE    Array[70, 100]
   VZ              DOUBLE    Array[70, 100]
   DENE            DOUBLE    Array[70, 100]
   ZEFF            DOUBLE    Array[70, 100]
```
