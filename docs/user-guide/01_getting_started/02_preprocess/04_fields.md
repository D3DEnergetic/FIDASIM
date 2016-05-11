title: Electromagnetic Fields

#Fields Structure

[TOC]

#Structure Variables
**time**: Time

* type: `Float64`
* rank: 0
* units: s

**data_source**: Source of electromagnetic fields data

* type: `String`
* rank: 0

**mask**: Boolean mask that indicates where the fields are well defined

* type: `Int16`
* rank: 2
* dims: [`nr`,`nz`]

**br**: Magnetic field in the r-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: T

**bt**: Magnetic field in the torodial-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: T

**bz**: Magnetic field in the z-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: T

**er**: Electric field in the r-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: V/m

**et**: Electric field in the torodial-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: V/m

**ez**: Electric field in the z-direction

* type: `Float64`
* rank: 2
* dims: [`nr`,`nz`]
* units: V/m

# Useful Routines

## extract_transp_geqdsk
`extract_transp_geqdsk` is a python routine that extracts a GEQDSK file out of a TRANSP run.
```
lstagner@computer:~/FIDASIM$ extract_transp_geqdsk -h
usage: extract_transp_geqdsk [-h] [-out OUTPUT_DIR] [-d DEVICE] [-y YEAR]
                             [-fid FILE_ID [FILE_ID ...]] [-t TIME [TIME ...]]
                             [-dt DELTA_TIME] [-bt BT_DIR] [-ip IP_DIR]
                             [-nr NUM_R] [-nz NUM_Z] [-nt NUM_THETA]
                             [-nv NUM_VERTICAL] [-nh NUM_HORIZONTAL]
                             [-nc NUM_CONTOURS] [-g GFILE] [-v] [-db]
                             directory runid

Extracts GEQDSK file from a TRANSP run using trxpl program

positional arguments:
  directory             Directory that contains TRANSP output files
  runid                 TRANSP run ID

optional arguments:
  -h, --help            show this help message and exit
  -out OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default = cwd
  -d DEVICE, --device DEVICE
                        Device Default = D3D
  -y YEAR, --year YEAR  Run year Default = Current Year
  -fid FILE_ID [FILE_ID ...], --file_id FILE_ID [FILE_ID ...]
                        File ID list e.g. 1 2 3 ...
  -t TIME [TIME ...], --time TIME [TIME ...]
                        Time of interest in seconds
  -dt DELTA_TIME, --delta_time DELTA_TIME
                        Averaging time window in seconds. Default = 0.005 s
  -bt BT_DIR, --bt_dir BT_DIR
                        Btorodial direction (1=ccw, -1=cw). Default = -1
  -ip IP_DIR, --ip_dir IP_DIR
                        Plasma current direction (1=ccw, -1=cw). Default = 1
  -nr NUM_R, --num_r NUM_R
                        Number of R points for cartesian overlay grid. Default
                        = 101
  -nz NUM_Z, --num_z NUM_Z
                        Number of Z points for cartesian overlay grid. Default
                        = 101
  -nt NUM_THETA, --num_theta NUM_THETA
                        Number of theta points for 2d splines. Default = 151
  -nv NUM_VERTICAL, --num_vertical NUM_VERTICAL
                        Number of vertical points for GEQDSK file. Default =
                        101
  -nh NUM_HORIZONTAL, --num_horizontal NUM_HORIZONTAL
                        Number of horizontal points for GEQDSK file. Default =
                        101
  -nc NUM_CONTOURS, --num_contours NUM_CONTOURS
                        Number of bdy contours for GEQDSK file. Default = 201
  -g GFILE, --gfile GFILE
                        GEQDSK file name. Default = "g[shot].[time]"
  -v, --verbose         Verbose
  -db, --debug          Debug mode
```

##read_geqdsk.pro
[read_geqdsk.pro](|url|/sourcefile/read_geqdsk.pro.html) is an IDL routine that reads a EFIT GEQDSK file and returns the fields structure.

##check_fields.pro
[check_fields.pro](|url|/sourcefile/check_fields.pro.html) is an IDL routine used internally by [PREFIDA](|url|/sourcefile/prefida.pro.html) to check if the fields structure has the correct format.

# Example Structure
```
IDL> help, fields
** Structure <456abb8>, 9 tags, length=350024, data length=350024, refs=1:
   TIME            DOUBLE           1.0000000
   DATA_SOURCE     STRING    '/u/lstagner/FIDASIM/test/g000001.01000'
   MASK            INT       Array[70, 100]
   BR              DOUBLE    Array[70, 100]
   BT              DOUBLE    Array[70, 100]
   BZ              DOUBLE    Array[70, 100]
   ER              DOUBLE    Array[70, 100]
   ET              DOUBLE    Array[70, 100]
   EZ              DOUBLE    Array[70, 100]
```
