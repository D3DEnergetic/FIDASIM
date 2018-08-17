title: Preprocessing Inputs

#Preprocessing Inputs

[TOC]

#Create FIDASIM input files using PREFIDA

FIDASIM requires inputs to be in a specific format.
PREFIDA([IDL](|url|/sourcefile/prefida.pro.html),[Python](|url|/sourcefile/preprocessing.py.html#prefida)) is an routine that takes the required inputs, checks their validity, and transforms them into a form FIDASIM understands.

PREFIDA is called as follows
```idl
IDL> prefida, inputs, grid, nbi, plasma, fields, fbm, spec=spec, npa=npa
```
or using Python
```python
>>> import preprocessing
>>> preprocessing.prefida(inputs, grid, nbi, plasma, fields, fbm, spec=spec, npa=npa)
```
where arguments are defined as follows. Click the argument's description for extreme detail.

* `inputs`: [General Settings](../03_technical/01_prefida_inputs.html#general-settings)
* `grid`: [Interpolation grid](../03_technical/01_prefida_inputs.html#interpolation-grid-structure)
* `nbi`: [Neutral Beam Geometry](../03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure)
* `fields`: [Electromagnetic Fields](../03_technical/01_prefida_inputs.html#fields-structure)
* `plasma`: [Plasma Parameters](../03_technical/01_prefida_inputs.html#plasma-structure)
* `fbm`: [Fast-ion Distribution](../03_technical/01_prefida_inputs.html#distribution-structure)
* `spec`: [Spectral Geometry](../03_technical/01_prefida_inputs.html#spectral-geometry-structure)
* `npa`: [NPA Geometry](../03_technical/01_prefida_inputs.html#npa-geometry-structure)

PREFIDA will create the following files

* Namelist File
* Geometry File
* Equilibrium File
* Distribution File

Most devices may have already setup helper routines to make running FIDASIM and Prefida easy. 
Click [here](./05_devices.html) to find out if someone has done your work for you.

#Making Grids
PREFIDA uses two types of grids: the [Interpolation Grid](../03_technical/01_prefida_inputs.html#interpolation-grid-structure) and the [Beam Grid](../03_technical/01_prefida_inputs.html#beam-grid-settings). 
The Interpolation Grid is a 3D cylindrical grid in R, Phi and Z. It is used for interpolating the [plasma parameters](../03_technical/01_prefida_inputs.html#plasma-structure) and the [electromagnetic fields](../03_technical/01_prefida_inputs.html#fields-structure). 
The routine `rz_grid`([IDL](|url|/sourcefile/rz_grid.pro.html),[Python](|url|/sourcefile/utils.py.html#rz_grid)) can be used to easily create the interpolation `grid` structure.

```idl
IDL> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz,phimin=phimin,phimax=phimax,nphi=nphi)
```
or in Python
```python
>>> from fidasim.utils import rz_grid
>>> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz,phimin=phimin,phimax=phimax,nphi=nphi)
```

The beam grid is a 3D grid used for most of the calculations in FIDASIM. It represents the 3D volume where the neutral beam lives and interacts with the plasma. 
To maximize the resolution of this grid it is useful to align the beam grid with the beam centerline.
The routine `beam_grid`([IDL](|url|/sourcefile/beam_grid.pro.html),[Python](|url|/sourcefile/utils.py.html#beam_grid)) calculates from the [neutral beam geometry](../03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure) the optimal [beam grid settings](../03_technical/01_prefida_inputs.html#beam-grid-settings) that would align the grid with the beam sightline.

```idl
IDL> bgrid = beam_grid(nbi,rstart)
```
or in Python
```python
>>> from fidasim.utils import beam_grid
>>> bgrid = beam_grid(nbi,rstart)
```

#Reading GEQDSK files
Most tokamaks use EFIT to reconstruct the MHD equilibrium.
To make things easy we provide the IDL routine [read_geqdsk.pro](|url|/sourcefile/read_geqdsk.pro.html) to calculate the [fields structure](../03_technical/01_prefida_inputs.html#fields-structure) from EFITs GEQDSK file.

```
IDL> fields = read_geqdsk('g159243.00300',grid,flux=flux)
```
where `grid` is the interpolation grid and the `flux` keyword is a named variable that recieves the torodial flux upon executation.

#Extracting GEQDSK file and Plasma Parameters from TRANSP
It is convenient to grab FIDASIM inputs from previously calculated TRANSP runs. 


The python script, `extract_transp_geqdsk`, can be used to extract the MHD equilibrium from TRANSP's `.DATA* files`.
For example:
```
extract_transp_geqdsk /p/transparch/result/NSTX/14 159243H06 
```
will create a GEQDSK file for every `.DATA*` file in the `159243H06` TRANSP run.
Run `extract_transp_geqdsk -h` for the full documentation.


The IDL routine [extract_transp_plasma.pro](|url|/sourcefile/extract_transp_plasma.pro.html) creates the [plasma structure](../03_technical/01_prefida_inputs.html#plasma-structure) at a given time. 

```
IDL> plasma = extract_transp_plasma("159243H06.CDF",1.02,grid,flux)
```
where `grid` is the interpolation grid and `flux` is the torodial flux.

#Translating NUBEAM Neutral Beam Geometry
The IDL routine [nubeam_geometry.pro](|url|/sourcefile/nubeam_geometry.pro.html) can be used to translate the NUBEAM neutral beam geometry definition into the [correct format](../03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure).
```
IDL> nbi = nubeam_geometry(nubeam)
```
where `nubeam` is a structure containing the NUBEAM geometry variables taken from the TRANSP namelist.

#Extracting the Fast-ion Distribution Function from TRANSP/NUBEAM
The python script, `extract_transp_fbm`, provides a easy way to extract the fast-ion distribution. For example: 
```
extract_transp_fbm /p/transparch/result/NSTX/14 159243H06
```
extracts a distribution function for every `.DATA*` file in the `159243H06` TRANSP run.
Run `extract_transp_fbm -h` for the full documentation.

#Reading NUBEAM/SPIRAL Fast-ion Distributions
Out of the box, FIDASIM provides IDL routines for reading different fast-ion distributions.
We provide routines for:

* [read_nubeam.pro](|url|/sourcefile/read_nubeam.pro.html): TRANSP/NUBEAM distribution functions
* [read_mc_nubeam.pro](|url|/sourcefile/read_mc_nubeam.pro.html): Monte Carlo TRANSP/NUBEAM Guiding Center distribution
* [read_spiral.pro](|url|/sourcefile/read_spiral.pro.html): SPIRAL Guiding Center distribution

```
IDL> f = read_nubeam(nubeam_distribution,grid,btipsign = -1) 
IDL> mcf = read_mc_nubeam(mc_nubeam_distribution,Ntotal=1e19,btipsign=-1)
IDL> s = read_spiral(spiral_file,Ntotal=1e19,btipsign=-1)
```
