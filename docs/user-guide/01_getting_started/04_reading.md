title: Reading Outputs

#Reading FIDASIM Outputs

[TOC]

##Using HDF5 tools
The HDF5 installation provides [several useful tools](https://www.hdfgroup.org/products/hdf5_tools/index.html#cmd) for inspecting HDF5 files.
```
lstagner@computer:~/TEST$ h5ls test_1a_geometry.h5
nbi                      Group
npa                      Group
spec                     Group
lstagner@computer:~/TEST$ h5ls test_1a_geometry.h5/nbi
axis                     Dataset {3}
data_source              Dataset {SCALAR}
divy                     Dataset {3}
divz                     Dataset {3}
focy                     Dataset {SCALAR}
focz                     Dataset {SCALAR}
name                     Dataset {SCALAR}
shape                    Dataset {SCALAR}
src                      Dataset {3}
widy                     Dataset {SCALAR}
widz                     Dataset {SCALAR}
lstagner@computer:~/TEST$ h5ls -v test_1a_geometry.h5/nbi/src
Opened "test_1a_geometry.h5" with sec2 driver.
src                      Dataset {3/3}
    Attribute: description scalar
        Type:      47-byte null-terminated ASCII string
        Data:  "Position of the center of the beam source grid"
    Attribute: units scalar
        Type:      3-byte null-terminated ASCII string
        Data:  "cm"
    Location:  1:5024
    Links:     1
    Chunks:    {3} 24 bytes
    Storage:   24 logical bytes, 19 allocated bytes, 126.32% utilization
    Filter-0:  shuffle-2 OPT {8}
    Filter-1:  deflate-1 OPT {9}
    Type:      native double
lstagner@computer:~/TEST$ h5ls -d test_1a_geometry.h5/nbi/src
src                      Dataset {3}
    Data:
        (0) 0, -502, 0
```

##Using IDL
The IDL routine [read_hdf5.pro](|url|/sourcefile/read_hdf5.pro.html) is provided to read any HDF5 file.

```idl
IDL> f = read_hdf5("./test_1a_geometry.h5")
IDL> help,f
** Structure <190eff78>, 4 tags, length=1848, data length=1830, refs=1:
   NBI             STRUCT    -> <Anonymous> Array[1]
   NPA             STRUCT    -> <Anonymous> Array[1]
   SPEC            STRUCT    -> <Anonymous> Array[1]
   DESCRIPTION     STRING    'Geometric quantities for FIDASIM'
IDL> help,f.nbi
** Structure <190bc848>, 13 tags, length=504, data length=498, refs=2:
   AXIS            STRUCT    -> <Anonymous> Array[1]
   DATA_SOURCE     STRUCT    -> <Anonymous> Array[1]
   DIVY            STRUCT    -> <Anonymous> Array[1]
   DIVZ            STRUCT    -> <Anonymous> Array[1]
   FOCY            STRUCT    -> <Anonymous> Array[1]
   FOCZ            STRUCT    -> <Anonymous> Array[1]
   NAME            STRUCT    -> <Anonymous> Array[1]
   SHAPE           STRUCT    -> <Anonymous> Array[1]
   SRC             STRUCT    -> <Anonymous> Array[1]
   WIDY            STRUCT    -> <Anonymous> Array[1]
   WIDZ            STRUCT    -> <Anonymous> Array[1]
   DESCRIPTION     STRING    'Neutral Beam Geometry'
   COORDINATE_SYSTEM
   STRING    'Right-handed cartesian'
IDL> help,f.nbi.src
** Structure <1909ade8>, 3 tags, length=56, data length=56, refs=2:
   DESCRIPTION     STRING    'Position of the center of the beam source grid'
   UNITS           STRING    'cm'
   DATA            DOUBLE    Array[3]
```

##Using Python
The [h5py](http://www.h5py.org) library can be used to read and write HDF5 files.

```python
>>> import h5py as h5
>>> f = h5.File("./test_1a_geometry.h5")
>>> [k for k in f.keys()]
['nbi', 'npa', 'spec']
>>> [k for k in f['nbi'].keys()]
['axis',
'data_source',
'divy',
'divz',
'focy',
'focz',
'name',
'shape',
'src',
'widy',
'widz']
>>> f["/nbi/src"].value
array([   0., -502.,    0.])
```

##Using Julia
The [HDF5.jl](https://github.com/JuliaLang/HDF5.jl) library can be used to read and write HDF5 files.
```julia
julia> using HDF5

julia> f = h5open("test_1a_geometry.h5")
HDF5 data file: test_1a_geometry.h5

julia> names(f)
3-element Array{ByteString,1}:
"nbi" 
"npa" 
"spec"

julia> names(f["/nbi"])
11-element Array{ByteString,1}:
"axis"       
"data_source"
"divy"       
"divz"       
"focy"       
"focz"       
"name"       
"shape"      
"src"        
"widy"       
"widz"       

julia> read(f["/nbi/src"])
3-element Array{Float64,1}:
0.0
-502.0
0.0

julia> h5read("test_1a_geometry.h5","/nbi/src")
3-element Array{Float64,1}:
0.0
-502.0
0.0
```
