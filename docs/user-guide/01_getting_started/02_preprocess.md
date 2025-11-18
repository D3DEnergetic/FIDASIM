title: Preprocessing Inputs

#Preprocessing Inputs

[TOC]

#Create FIDASIM input files using PREFIDA

FIDASIM requires inputs to be in a specific format.
PREFIDA([IDL](|url|/sourcefile/prefida.pro.html),[Python](|url|/sourcefile/preprocessing.py.html#prefida)) is an routine that takes the required inputs, checks their validity, and transforms them into a form FIDASIM understands.

PREFIDA is called as follows
```idl
IDL> prefida, inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa
```
or using Python
```python
>>> import preprocessing
>>> preprocessing.prefida(inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa)
```
where arguments are defined as follows. Click the argument description for extreme detail.

* `inputs`: [General Settings](../03_technical/01_prefida_inputs.html#general-settings)
* `grid`: [Interpolation grid](../03_technical/01_prefida_inputs.html#interpolation-grid-structure)
* `nbi`: [Neutral Beam Geometry](../03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure)
* `fields`: [Electromagnetic Fields](../03_technical/01_prefida_inputs.html#fields-structure)
* `plasma`: [Plasma Parameters](../03_technical/01_prefida_inputs.html#plasma-structure)
* `dist`: [Fast-ion Distribution](../03_technical/01_prefida_inputs.html#distribution-structure)
* `spec`: [Spectral Geometry](../03_technical/01_prefida_inputs.html#spectral-geometry-structure)
* `npa`: [NPA Geometry](../03_technical/01_prefida_inputs.html#npa-geometry-structure)

## Including Molecular Hydrogen (H2) for Passive Signals

To include molecular hydrogen in passive signal calculations, add the `denm` field to your plasma structure:

```python
# Example: Adding molecular hydrogen density
plasma['denm'] = denm_array  # Shape: [nthermal, nr, nz] or [nthermal, nr, nz, nphi]
```

Where `denm` is the molecular hydrogen density in cm^-3. This is optional - if not provided, only atomic charge exchange will be considered.

PREFIDA will create the following files

* Namelist File
* Geometry File
* Equilibrium File
* Distribution File

Most devices may have already setup helper routines to make running FIDASIM and Prefida easy. 
Click [here](./05_devices.html) to find out if someone has done your work for you.

#Making Grids

PREFIDA uses two types of grids**: the [Interpolation Grid](../03_technical/01_prefida_inputs.html#interpolation-grid-structure) and the [Beam Grid](../03_technical/01_prefida_inputs.html#beam-grid-settings). 
By default, axisymmetry is assumed and the Interpolation Grid is a 2D grid in the R-Z plane that is used for interpolating the [plasma parameters](../03_technical/01_prefida_inputs.html#plasma-structure) and the [electromagnetic fields](../03_technical/01_prefida_inputs.html#fields-structure).
A 3D cylindrical grid in R, Z and Phi can be created if the user inputs phi variable information. 
The routine `rz_grid`([IDL](|url|/sourcefile/rz_grid.pro.html),[Python](|url|/sourcefile/utils.py.html#rz_grid)) accomplishes the task and creates the `grid` structure. For example, the command below will create a 2D grid,
```idl
IDL> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz)
```
whereas the following command will create a 3D grid,
```idl
IDL> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz,phimin=phimin,phimax=phimax,nphi=nphi)
```
In Python, the 2D grid can be created with,
```python
>>> from fidasim.utils import rz_grid
>>> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz)
```
and the 3D grid with,
```python
>>> from fidasim.utils import rz_grid
>>> grid = rz_grid(rmin,rmax,nr,zmin,zmax,nz,phimin=phimin,phimax=phimax,nphi=nphi)
```
The output 2D grid structure will have Phi = 0.0 and nphi = 1, but the 3D grid structure will have values based on what the user input.

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

FIDASIM creates a third grid for [passive calculations](../02_physics/08_passive.html#passive-grid).

#Reading GEQDSK files
Most tokamaks use EFIT to reconstruct the MHD equilibrium.
To make things easy we provide the IDL routine [read_geqdsk.pro](|url|/sourcefile/read_geqdsk.pro.html) to calculate the [fields structure](../03_technical/01_prefida_inputs.html#fields-structure) from EFITs GEQDSK file.

```idl
IDL> fields = read_geqdsk('g159243.00300',grid,rho=rho,btipsign=btipsign)
```
or in Python
```python
>>> from fidasim.utils import read_geqdsk
>>> fields, rho, btipsign = read_geqdsk('g159243.00300',grid)
```
where `grid` is the interpolation grid, `rho` keyword is a named variable that recieves the sqrt(normalized torodial flux) upon executation, and `btipsign` is a named variable that recieves the Bt-Ip sign (-1 for anti-parallel, 1 for parallel).

#Extracting GEQDSK file and Plasma Parameters from TRANSP
It is convenient to grab FIDASIM inputs from previously calculated TRANSP runs. 


The python script, `extract_transp_geqdsk`, can be used to extract the MHD equilibrium from TRANSP's `.DATA* files`.
For example:
```bash
extract_transp_geqdsk /p/transparch/result/NSTX/14 159243H06 
```
will create a GEQDSK file for every `.DATA*` file in the `159243H06` TRANSP run.
Run `extract_transp_geqdsk -h` for the full documentation.


The IDL routine [extract_transp_plasma.pro](|url|/sourcefile/extract_transp_plasma.pro.html) or the equivalent Python function creates the [plasma structure](../03_technical/01_prefida_inputs.html#plasma-structure) at a given time. 

```idl
IDL> plasma = extract_transp_plasma("159243H06.CDF",1.02,grid,rho)
```
or in Python
```python
>>> from fidasim.utils import extract_transp_plasma
>>> plasma = extract_transp_plasma("159243H06.CDF",1.02,grid,rho)
```
where `grid` is the interpolation grid and `rho` is the sqrt(normalized torodial flux).

#Translating NUBEAM Neutral Beam Geometry
The IDL routine [nubeam_geometry.pro](|url|/sourcefile/nubeam_geometry.pro.html) can be used to translate the NUBEAM neutral beam geometry definition into the [correct format](../03_technical/01_prefida_inputs.html#neutral-beam-geometry-structure).
```idl
IDL> nbi = nubeam_geometry(nubeam)
```
or in Python
```python
>>> from fidasim.utils import nubeam_geometry
>>> nbi = nubeam_geometry(nubeam)
```
where `nubeam` is a structure/dictionary containing the NUBEAM geometry variables taken from the TRANSP namelist.

#Extracting the Fast-ion Distribution Function from TRANSP/NUBEAM
The python script, `extract_transp_fbm`, provides a easy way to extract the fast-ion distribution. For example:
```bash
extract_transp_fbm /p/transparch/result/NSTX/14 159243H06
```
extracts a distribution function for every `.DATA*` file in the `159243H06` TRANSP run.
Run `extract_transp_fbm -h` for the full documentation.

#Reading NUBEAM/SPIRAL Fast-ion Distributions
Out of the box, FIDASIM provides IDL and Python routines for reading different fast-ion distributions.
We provide routines for:

* [read_nubeam.pro](|url|/sourcefile/read_nubeam.pro.html): TRANSP/NUBEAM distribution functions
* [read_mc_nubeam.pro](|url|/sourcefile/read_mc_nubeam.pro.html): Monte Carlo TRANSP/NUBEAM Guiding Center distribution
* [read_spiral.pro](|url|/sourcefile/read_spiral.pro.html): SPIRAL Guiding Center distribution

```idl
IDL> f = read_nubeam(nubeam_distribution_file,grid,btipsign = -1) 
IDL> mcf = read_mc_nubeam(mc_nubeam_distribution_file,Ntotal=1e19,btipsign=-1)
IDL> s = read_spiral(spiral_file,Ntotal=1e19,btipsign=-1)
```
or in Python
```python
>>> from fidasim.utils import read_nubeam
>>> f = read_nubeam(nubeam_distribution_file, grid, btipsign=-1)
```
