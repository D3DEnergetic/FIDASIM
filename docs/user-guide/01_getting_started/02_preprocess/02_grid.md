title: Interpolation Grid

#Grid Structure
The `grid` structure contains the definition of the 2D R-Z grid that the [plasma parameters](./05_plasma.html) and [electromagnetic fields](./04_fields.html) are mapped onto. 

[TOC]

#Useful Routines
##rz_grid.pro
An IDL routine [rz_grid.pro](|url|/sourcefile/rz_grid.pro.html) is provided to easily define the interpolation grid.

##check_grid.pro
[check_grid.pro](|url|/sourcefile/check_grid.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the grid structure has the correct format.

#Example Structure
```
IDL> help, grid
** Structure <44c8108>, 6 tags, length=113368, data length=113364, refs=1:
   R2D             DOUBLE    Array[70, 100]
   Z2D             DOUBLE    Array[70, 100]
   R               DOUBLE    Array[70]
   Z               DOUBLE    Array[100]
   NR              INT             70
   NZ              INT            100
```

#Structure Variables
**nr**: Number of radii

* type: `Int16`
* rank: 0

**nz**: Number of z values

* type: `Int16`
* rank: 0

**r2d**: 2D array of radii

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm

**z2d**: 2D array of z values [cm]

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: cm

**r**: Array of radii [cm]

* type: `Float64`
* rank: 1
* dims: [`nr`]
* units: cm

**z**: Array of z values [cm]

* type: `Float64`
* rank: 1
* dims: [`nz`]
* units: cm

