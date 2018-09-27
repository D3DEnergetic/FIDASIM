title: Reading Outputs

#Reading FIDASIM Outputs

[TOC]

##Using HDF5 tools
The HDF5 installation provides [several useful tools](https://www.hdfgroup.org/products/hdf5_tools/index.html#cmd) for inspecting HDF5 files.
```
lstagner@computer:~/TEST$ h5ls test_geometry.h5
nbi                      Group
npa                      Group
spec                     Group
lstagner@computer:~/TEST$ h5ls test_geometry.h5/nbi
adist                    Dataset {1}
aoffy                    Dataset {1}
aoffz                    Dataset {1}
ashape                   Dataset {1}
awidy                    Dataset {1}
awidz                    Dataset {1}
axis                     Dataset {3}
data_source              Dataset {SCALAR}
divy                     Dataset {3}
divz                     Dataset {3}
focy                     Dataset {SCALAR}
focz                     Dataset {SCALAR}
name                     Dataset {SCALAR}
naperture                Dataset {SCALAR}
shape                    Dataset {SCALAR}
src                      Dataset {3}
widy                     Dataset {SCALAR}
widz                     Dataset {SCALAR}
lstagner@computer:~/TEST$ h5ls -v test_geometry.h5/nbi/src
Opened "test_geometry.h5" with sec2 driver.
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
    Storage:   24 logical bytes, 20 allocated bytes, 120.00% utilization
    Filter-0:  shuffle-2 OPT {8}
    Filter-1:  deflate-1 OPT {9}
    Type:      native double
lstagner@computer:~/TEST$ h5ls -d test_geometry.h5/nbi/src
src                      Dataset {3}
    Data:
        (0) 0, -532, 0
```

##Using IDL
The IDL routine [read_hdf5.pro](|url|/sourcefile/read_hdf5.pro.html) is provided to read any HDF5 file.

```idl
IDL> f = read_hdf5("./test_geometry.h5")
IDL> help,f
** Structure <284b698>, 4 tags, length=2224, data length=2194, refs=1:
   NBI             STRUCT    -> <Anonymous> Array[1]
   NPA             STRUCT    -> <Anonymous> Array[1]
   SPEC            STRUCT    -> <Anonymous> Array[1]
   DESCRIPTION     STRING    'Geometric quantities for FIDASIM'
IDL> help,f.nbi
** Structure <27fb4d8>, 20 tags, length=752, data length=734, refs=2:
   ADIST           STRUCT    -> <Anonymous> Array[1]
   AOFFY           STRUCT    -> <Anonymous> Array[1]
   AOFFZ           STRUCT    -> <Anonymous> Array[1]
   ASHAPE          STRUCT    -> <Anonymous> Array[1]
   AWIDY           STRUCT    -> <Anonymous> Array[1]
   AWIDZ           STRUCT    -> <Anonymous> Array[1]
   AXIS            STRUCT    -> <Anonymous> Array[1]
   DATA_SOURCE     STRUCT    -> <Anonymous> Array[1]
   DIVY            STRUCT    -> <Anonymous> Array[1]
   DIVZ            STRUCT    -> <Anonymous> Array[1]
   FOCY            STRUCT    -> <Anonymous> Array[1]
   FOCZ            STRUCT    -> <Anonymous> Array[1]
   NAME            STRUCT    -> <Anonymous> Array[1]
   NAPERTURE       STRUCT    -> <Anonymous> Array[1]
   SHAPE           STRUCT    -> <Anonymous> Array[1]
   SRC             STRUCT    -> <Anonymous> Array[1]
   WIDY            STRUCT    -> <Anonymous> Array[1]
   WIDZ            STRUCT    -> <Anonymous> Array[1]
   DESCRIPTION     STRING    'Neutral Beam Geometry'
   COORDINATE_SYSTEM
                   STRING    'Right-handed cartesian'
IDL> help,f.nbi.src
** Structure <27f5908>, 3 tags, length=56, data length=56, refs=2:
   DESCRIPTION     STRING    'Position of the center of the beam source grid'
   UNITS           STRING    'cm'
   DATA            DOUBLE    Array[3]
```

##Using Python
The [h5py](http://www.h5py.org) library can be used to read and write HDF5 files.

```python
>>> import h5py as h5
>>> f = h5.File("./test_geometry.h5")
>>> [k for k in f.keys()]
['nbi', 'npa', 'spec']
>>> [k for k in f['nbi'].keys()]
['adist',
 'aoffy',
 'aoffz',
 'ashape',
 'awidy',
 'awidz',
 'axis',
 'data_source',
 'divy',
 'divz',
 'focy',
 'focz',
 'name',
 'naperture',
 'shape',
 'src',
 'widy',
 'widz']
>>> f["/nbi/src"].value
array([   0., -532.,    0.])
```

##Using Julia
The [HDF5.jl](https://github.com/JuliaLang/HDF5.jl) library can be used to read and write HDF5 files.
```julia
julia> using HDF5

julia> f = h5open("test_geometry.h5")
HDF5 data file: test_geometry.h5

julia> names(f)
3-element Array{ByteString,1}:
"nbi" 
"npa" 
"spec"

julia> names(f["/nbi"])
11-element Array{ByteString,1}:
"adist"
"aoffy"
"aoffz"
"ashape"
"awidy"
"awidz"
"axis"
"data_source"
"divy"
"divz"
"focy"
"focz"
"name"
"naperture"
"shape"
"src"
"widy"
"widz"

julia> read(f["/nbi/src"])
3-element Array{Float64,1}:
0.0
-532.0
0.0

julia> h5read("test_geometry.h5","/nbi/src")
3-element Array{Float64,1}:
0.0
-532.0
0.0
```

#Visualization: Outputs
Visualizing your outputs can be done by executing `plot_outputs` found in `lib/scripts/`

You can plot your spectra and npa files, and print the neutron rate in various ways.
Let's start with the easiest method first.
Assuming all of your output files are located in the same folder, the following command will search that folder for all run IDs present.
Then, it will plot all of the spectra and npa data for each run ID on every FIDA and NPA channel.
It will also print out the neutron rate.
```
plot_outputs -d /p/fida/lstagner/TEST/ -s -as -n -an
```
If you are interested in only one shot, add `-r` to your command to indicate the run ID you wish to analyze.

What if your files are scattered around in different folders?
What if you are only interested in the active FIDA emission?
What if you want to view channels 1 and 3?
Don't worry, this is also easy to accomplish.
The example below does this for two different files located in different folders.
```
plot_outputs -p /p/fida/lstagner/TEST/test_spectra.h5 /p/fida/lstagner/DIFFERENT_TEST/different_test_spectra.h5 -f -ls 1 3
```

Now that you have a grasp on how the script works, feel free to look at the help documentation to see what else the code is capable of.
