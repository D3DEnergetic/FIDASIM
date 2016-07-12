title: Plasma Parameters, Fields, and Distributions

[TOC]

---

<style>
table {
width: 100%;
}
table,th,td {
border: 1px solid black;
border-collapse: collapse;
}
th, td {
padding: 5px;
}
th {
text-align: center;
}
</style>

# Plasma Parameters and Fields
![](|media|/cross_section.png){: style="float: left"}
It is common in Tokamaks to input the plasma parameters and fields as flux functions.
However, this approach makes several assumptions that may not be true in non-tokamak devices.
In FIDASIM the plasma parameters and fields are mapped onto a 2D R-Z grid.
This only assumes axisymmetric profiles within the region of interest.
In addition to plasma parameters a boolean mask is used to indicate where the plasma and fields are well defined.


FIDASIM reads the plasma parameters and fields from an HDF5 file structured as follows...
```
[runid]_equilibrium.h5
├── plasma          
└── fields
```
where the `plasma` group has the following datasets

|       Variable      |   Type  | Rank |  Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:-----------------------------------------|
| `nr`                | Int16   | 0    | NA          | NA    | Number of radii                          |
| `nz`                | Int16   | 0    | NA          | NA    | Number of z values                       |
| `r`                 | Float64 | 1    | [`nr`]      | cm    | Array of radii                           |
| `z`                 | Float64 | 1    | [`nz`]      | cm    | Array of z values                        |
| `r2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of radii `r = r2d(r,z)`         |
| `z2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of z values `z = z2d(r,z)`      |
| `time`              | Float64 | 0    | NA          | s     | Time when the plasma parameter data was collected             |
| `data_source`       | String  | 0    | NA          | NA    | Source of the plasma parameter data                           |
| `mask`              | Int16   | 2    | [`nr`,`nz`] | NA    | Boolean mask that indicates where the plasma is well defined  |
| `te`                | Float64 | 2    | [`nr`,`nz`] | keV   | Electron temperature                                          |
| `ti`                | Float64 | 2    | [`nr`,`nz`] | keV   | Ion temperature                                               |
| `dene`              | Float64 | 2    | [`nr`,`nz`] | cm^-3 | Electron density                                              |
| `zeff`              | Float64 | 2    | [`nr`,`nz`] | NA    | Z-effective                                                   |
| `vr`                | Float64 | 2    | [`nr`,`nz`] | cm/s  | Radial component of the bulk plasma rotation/flow             |
| `vt`                | Float64 | 2    | [`nr`,`nz`] | cm/s  | Torodial/Phi component of the bulk plasma rotation/flow       |
| `vz`                | Float64 | 2    | [`nr`,`nz`] | cm/s  | Z component of the bulk plasma rotation/flow                  |

and where the `fields` group has the following datasets

|       Variable      |   Type  | Rank |  Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:-----------------------------------------|
| `nr`                | Int16   | 0    | NA          | NA    | Number of radii                          |
| `nz`                | Int16   | 0    | NA          | NA    | Number of z values                       |
| `r`                 | Float64 | 1    | [`nr`]      | cm    | Array of radii                           |
| `z`                 | Float64 | 1    | [`nz`]      | cm    | Array of z values                        |
| `r2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of radii `r = r2d(r,z)`         |
| `z2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of z values `z = z2d(r,z)`      |
| `time`              | Float64 | 0    | NA          | s     | Time when the fields data were collected/reconstructed        |
| `data_source`       | String  | 0    | NA          | NA    | Source of the fields data                                     |
| `mask`              | Int16   | 2    | [`nr`,`nz`] | NA    | Boolean mask that indicates where the fields are well defined |
| `br`                | Float64 | 2    | [`nr`,`nz`] | T     | Radial component of the magnetic field                        |
| `bt`                | Float64 | 2    | [`nr`,`nz`] | T     | Torodial/Phi component of the magnetic field                  |
| `bz`                | Float64 | 2    | [`nr`,`nz`] | T     | Z component of the magnetic field                             |
| `er`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Radial component of the electric field                        |
| `et`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Torodial/Phi component of the electric field                  |
| `ez`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Z component of the electric field                             |

# Distributions
The fast-ion distribution function typically does not have a functional form.
This means that we need to input the fast-ion distribution function into FIDASIM.
FIDASIM accepts three different types of fast-ion distributions.

1. Guiding Center Distribution Function: F(E,p,R,Z)
2. Guiding Center Monte Carlo distribution
3. Full-orbit Monte Carlo distribution

Supporting different types of distributions greatly facilitates the study of fast-ion transport.

FIDASIM reads the fast-ion distribution from an HDF5 file that has the following variables depending on the distribution type.

## Guiding Center Distribution Function: F(E,p,R,Z)

|       Variable      |   Type  | Rank |           Dimensions           |          Units         |           Description           |
|:-------------------:|:-------:|:----:|:------------------------------:|:----------------------:|:--------------------------------| 
| `type`              | Int16   | 0    | NA                             | NA                     | Distribution type (1)           |
| `nr`                | Int16   | 0    | NA                             | NA                     | Number of radii                 |
| `nz`                | Int16   | 0    | NA                             | NA                     | Number of z values              |
| `r`                 | Float64 | 1    | [`nr`]                         | cm                     | Array of radii                  |
| `z`                 | Float64 | 1    | [`nz`]                         | cm                     | Array of z values               |
| `r2d`               | Float64 | 2    | [`nr`,`nz`]                    | cm                     | 2D array of radii `r = r2d(r,z)`|
| `z2d`               | Float64 | 2    | [`nr`,`nz`]                    | cm                     | 2D array of z values `z = z2d(r,z)`|
| `time`              | Float64 | 0    | NA                             | s                      | Time of the distribution        |
| `data_source`       | String  | 0    | NA                             | NA                     | Source of the distribution data |
| `nenergy`           | Int16   | 0    | NA                             | NA                     | Number of energy values         |
| `npitch`            | Int16   | 0    | NA                             | NA                     | Number of pitch values          |
| `energy`            | Float64 | 1    | [`nenergy`]                    | keV                    | Energy array                    |
| `pitch`             | Float64 | 1    | [`npitch`]                     | NA                     | Pitch array w.r.t magnetic field|
| `denf`              | Float64 | 2    | [`nr`,`nz`]                    | cm^-3                  | Fast-ion density                |
| `f`                 | Float64 | 4    | [`nenergy`,`npitch`,`nr`,`nz`] | fast-ions/(dE dP cm^3) | Fast-ion distribution F(E,p,R,Z)|

The Guiding Center Distribution Function is a function of Energy, pitch, R, and Z.
The distribution is mapped on the same R-Z grid as the plasma parameters and fields.

## Guiding Center Monte Carlo Distribution

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `type`              | Int16   | 0    | NA           | NA    | Distribution type (2)           |
| `time`              | Float64 | 0    | NA           | s     | Time of the distribution        |
| `data_source`       | String  | 0    | NA           | NA    | Source of the distribution data |
| `nparticle`         | Int32   | 0    | NA           | NA    | Number of MC particles          |
| `nclass`            | Int16   | 0    | NA           | NA    | Number of orbit classes         |
| `class`             | Int16   | 1    | [`nparticle`]| NA    | Orbit class of the MC particle  |
| `weight`            | Float64 | 1    | [`nparticle`]| fast-ions| Weight of the MC particle    |
| `r`                 | Float64 | 1    | [`nparticle`]| cm    | R positions of the MC particle  |
| `z`                 | Float64 | 1    | [`nparticle`]| cm    | Z positions of the MC particle  |
| `energy`            | Float64 | 1    | [`nparticle`]| keV   | Energy of the MC particle       |
| `pitch`             | Float64 | 1    | [`nparticle`]| NA    | Pitch w.r.t the magnetic field of the MC particle |

The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles.

The `class` variable can take values in the range of 1:`nclass`. 
If there are multiple classes of particles the FIDA signal for each class will be calculated.

## Full-orbit Monte Carlo Distribution

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `type`              | Int16   | 0    | NA           | NA    | Distribution type (3)           |
| `time`              | Float64 | 0    | NA           | s     | Time of the distribution        |
| `data_source`       | String  | 0    | NA           | NA    | Source of the distribution data |
| `nparticle`         | Int32   | 0    | NA           | NA    | Number of MC particles          |
| `nclass`            | Int16   | 0    | NA           | NA    | Number of orbit classes         |
| `class`             | Int16   | 1    | [`nparticle`]| NA    | Orbit class of the MC particle  |
| `weight`            | Float64 | 1    | [`nparticle`]| fast-ions| Weight of the MC particle    |
| `r`                 | Float64 | 1    | [`nparticle`]| cm    | R positions of the MC particle  |
| `z`                 | Float64 | 1    | [`nparticle`]| cm    | Z positions of the MC particle  |
| `vr`                | Float64 | 1    | [`nparticle`]| cm/s  | Radial component of the MC particle velocity |
| `vt`                | Float64 | 1    | [`nparticle`]| cm/s  | Torodial/Phi component of the MC particle velocity |
| `vz`                | Float64 | 1    | [`nparticle`]| cm/s  | Z component of the MC particle velocity |

The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles.

The `class` variable can take values in the range of 1:`nclass`.
If there are multiple classes of particles the FIDA signal for each class will be calculated.

# Relevent Namelist Settings
* `ai`: Ion mass [amu]
* `impurity_charge`: Impurity charge number 5=Boron, 6=Carbon, ...
* `ab`: Fast/Beam-ion mass [amu]
* `equilibrium_file`: Equilibrium file location
* `distribution_file`: Distribution file location

# Fortran References

* [[InterpolationGrid]]: Definition of R-Z grid
* [[Profiles]] and [[LocalProfiles]]: Derived type for Plasma parameters
* [[EMFields]] and [[LocalEMFields]]: Derived type for Fields
* [[get_plasma]]: Gets plasma parameters at a given position
* [[get_fields]]: Gets fields at a given position
* [[get_distribution]]: Gets fast-ion distribution at a given position
* [[read_equilibrium]]: Reads equilbrium file into [[Equilibrium]] structure
* [[FastIonDistribution]]: Derived type for describing GC distribution functions
* [[FastIon]]: Derived type for describing a Monte Carlo particle
* [[FastIonParticles]]: Defines a distribution of Monte Carlo particles
* [[read_distribution]]: Reads distribution file
