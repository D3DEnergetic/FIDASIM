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
The default settings of FIDASIM are to map the plasma parameters and fields onto the 2D R-Z grid. If the user inputs any Phi variable information, the plasma parameters and fields will be mapped onto a 3D R-Z-Phi cylindrical grid. Furthermore, a boolean mask is used to indicate where the plasma and fields are well defined.


FIDASIM reads the plasma parameters and fields from an HDF5 file structured as follows...
```
[runid]_equilibrium.h5
├── plasma
└── fields
```
where the `plasma` group has the following datasets

|       Variable      |   Type  | Rank |            Dimensions           | Units |               Description                |
|:-------------------:|:-------:|:----:|:-------------------------------:|:-----:|:-----------------------------------------|
| `nr`                | Int16   | 0    | NA                              | NA    | Number of radii                          |
| `nz`                | Int16   | 0    | NA                              | NA    | Number of z values                       |
| `nphi`              | Int16   | 0    | NA                              | NA    | Number of phi values (Optional)          |
| `nthermal`          | Int16   | 0    | NA                              | NA    | Number of hydrogenic, thermal ion species|
| `impurity_charge`   | Int16   | 0    | NA                              | NA    | Main impurity charge number              |
| `species_mass`      | Float64 | 1    | [`nthermal`]                    | amu   | Thermal ion species mass                 |
| `r`                 | Float64 | 1    | [`nr`]                          | cm    | Array of radii                           |
| `z`                 | Float64 | 1    | [`nz`]                          | cm    | Array of z values                        |
| `phi`               | Float64 | 1    | [`nphi`]                        | rad   | Array of phi values (Optional)           |
| `r2d`               | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm    | 2D/3D array of radii `r = r2d(r,z[,phi])`|
| `z2d`               | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm    | 2D/3D array of z values `z = z2d(r,z,[phi])`|
| `time`              | Float64 | 0    | NA                              | s     | Time when the plasma parameter data was collected            |
| `data_source`       | String  | 0    | NA                              | NA    | Source of the plasma parameter data                          |
| `mask`              | Int16   | 2/3  | [`nr`,`nz`[,`nphi`]]            | NA    | Boolean mask that indicates where the plasma is well defined |
| `te`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | keV   | Electron temperature                                         |
| `ti`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | keV   | Ion temperature                                              |
| `dene`              | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm^-3 | Electron density                                             |
| `denimp`            | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm^-3 | Impurity density                                             |
| `deni`              | Float64 | 3/4  | [`nthermal`,`nr`,`nz`[,`nphi`]] | cm^-3 | Ion density for each thermal species                         |
| `denn`              | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm^-3 | Cold neutral density                                         |
| `zeff`              | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | NA    | Z-effective                                                  |
| `vr`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm/s  | Radial component of the bulk plasma rotation/flow            |
| `vt`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm/s  | Torodial/Phi component of the bulk plasma rotation/flow      |
| `vz`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]]            | cm/s  | Z component of the bulk plasma rotation/flow                 |
| `description`       | String  | 0    | NA                              | NA    | Plasma Parameters                        |
| `coordinate system` | String  | 0    | NA                              | NA    | Cylindrical                              |

and where the `fields` group has the following datasets

|       Variable      |   Type  | Rank |       Dimensions     | Units |               Description                   |
|:-------------------:|:-------:|:----:|:--------------------:|:-----:|:--------------------------------------------|
| `nr`                | Int16   | 0    | NA                   | NA    | Number of radii                             |
| `nz`                | Int16   | 0    | NA                   | NA    | Number of z values                          |
| `nphi`              | Int16   | 0    | NA                   | NA    | Number of phi values (Optional)             |
| `r`                 | Float64 | 1    | [`nr`]               | cm    | Array of radii                              |
| `z`                 | Float64 | 1    | [`nz`]               | cm    | Array of z values                           |
| `phi`               | Float64 | 1    | [`nphi`]             | rad   | Array of phi values (Optional)              |
| `r2d`               | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | cm    | 2D/3D array of radii `r = r2d(r,z[,phi])`   |
| `z2d`               | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | cm    | 2D/3D array of z values `z = z2d(r,z,[phi])`|
| `time`              | Float64 | 0    | NA                   | s     | Time when the fields data were collected/reconstructed        |
| `data_source`       | String  | 0    | NA                   | NA    | Source of the fields data                                     |
| `mask`              | Int16   | 2/3  | [`nr`,`nz`[,`nphi`]] | NA    | Boolean mask that indicates where the fields are well defined |
| `br`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | T     | Radial component of the magnetic field                        |
| `bt`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | T     | Torodial/Phi component of the magnetic field                  |
| `bz`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | T     | Z component of the magnetic field                             |
| `er`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | V/m   | Radial component of the electric field                        |
| `et`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | V/m   | Torodial/Phi component of the electric field                  |
| `ez`                | Float64 | 2/3  | [`nr`,`nz`[,`nphi`]] | V/m   | Z component of the electric field                             |
| `description`       | String  | 0    | NA                   | NA    | Electromagnetic Field                       |
| `coordinate system` | String  | 0    | NA                   | NA    | Cylindrical                                 |

# Distributions
The fast-ion distribution function typically does not have a functional form.
This means that we need to input the fast-ion distribution function into FIDASIM.
FIDASIM accepts three different types of fast-ion distributions.

1. Guiding Center Distribution Function: F(E,p,R,Z[,Phi])
2. Guiding Center Monte Carlo distribution
3. Full-orbit Monte Carlo distribution

Supporting different types of distributions greatly facilitates the study of fast-ion transport.

FIDASIM reads the fast-ion distribution from an HDF5 file that has the following variables depending on the distribution type.

## Guiding Center Distribution Function: F(E,p,R,Z[,Phi])

|       Variable      |   Type  | Rank |           Dimensions           |          Units         |           Description           |
|:-------------------:|:-------:|:----:|:------------------------------:|:----------------------:|:--------------------------------|
| `type`              | Int16   | 0    | NA                                      | NA                     | Distribution type (1)                  |
| `nr`                | Int16   | 0    | NA                                      | NA                     | Number of radii                        |
| `nz`                | Int16   | 0    | NA                                      | NA                     | Number of z values                     |
| `nphi`              | Int16   | 0    | NA                                      | NA                     | Number of phi values (Optional)        |
| `r`                 | Float64 | 1    | [`nr`]                                  | cm                     | Array of radii                         |
| `z`                 | Float64 | 1    | [`nz`]                                  | cm                     | Array of z values                      |
| `phi`               | Float64 | 1    | [`nphi`]                                | cm                     | Array of phi values (Optional)         |
| `r2d`               | Float64 | 2    | [`nr`,`nz`]                             | cm                     | 2D array of radii `r = r2d(r,z)`       |
| `z2d`               | Float64 | 2    | [`nr`,`nz`]                             | cm                     | 2D array of z values `z = z2d(r,z)`    |
| `time`              | Float64 | 0    | NA                                      | s                      | Time of the distribution               |
| `data_source`       | String  | 0    | NA                                      | NA                     | Source of the distribution data        |
| `nenergy`           | Int16   | 0    | NA                                      | NA                     | Number of energy values                |
| `npitch`            | Int16   | 0    | NA                                      | NA                     | Number of pitch values                 |
| `energy`            | Float64 | 1    | [`nenergy`]                             | keV                    | Energy array                           |
| `pitch`             | Float64 | 1    | [`npitch`]                              | NA                     | Pitch array w.r.t magnetic field       |
| `denf`              | Float64 | 3    | [`nr`,`nz`[,`nphi`]]                      | cm^-3                  | Fast-ion density                       |
| `f`                 | Float64 | 5    | [`nenergy`,`npitch`,`nr`,`nz`[,`nphi`]] | fast-ions/(dE dP cm^3) | Fast-ion distribution F(E,p,R,Z[,Phi]) |

The Guiding Center Distribution Function is a function of Energy, pitch, R, Z and Phi.
The distribution can be mapped onto the 2D R-Z grid or 3D cylindrical grid, where the plasma parameters and fields are defined.

## Guiding Center Monte Carlo Distribution

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------|
| `type`              | Int16   | 0    | NA           | NA    | Distribution type (2)                             |
| `time`              | Float64 | 0    | NA           | s     | Time of the distribution                          |
| `data_source`       | String  | 0    | NA           | NA    | Source of the distribution data                   |
| `nparticle`         | Int32   | 0    | NA           | NA    | Number of MC particles                            |
| `nclass`            | Int16   | 0    | NA           | NA    | Number of orbit classes                           |
| `class`             | Int16   | 1    | [`nparticle`]| NA    | Orbit class of the MC particle                    |
| `weight`            | Float64 | 1    | [`nparticle`]| fast-ions| Weight of the MC particle                      |
| `r`                 | Float64 | 1    | [`nparticle`]| cm    | R positions of the MC particle                    |
| `z`                 | Float64 | 1    | [`nparticle`]| cm    | Z positions of the MC particle                    |
| `phi`               | Float64 | 1    | [`nparticle`]| rad   | Phi positions of the MC particle (Optional)       |
| `energy`            | Float64 | 1    | [`nparticle`]| keV   | Energy of the MC particle                         |
| `pitch`             | Float64 | 1    | [`nparticle`]| NA    | Pitch w.r.t the magnetic field of the MC particle |

The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles.

The `class` variable can take values in the range of 1:`nclass`.
If there are multiple classes of particles the FIDA signal for each class will be calculated.

## Full-orbit Monte Carlo Distribution

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------|
| `type`              | Int16   | 0    | NA           | NA    | Distribution type (3)                       |
| `time`              | Float64 | 0    | NA           | s     | Time of the distribution                    |
| `data_source`       | String  | 0    | NA           | NA    | Source of the distribution data             |
| `nparticle`         | Int32   | 0    | NA           | NA    | Number of MC particles                      |
| `nclass`            | Int16   | 0    | NA           | NA    | Number of orbit classes                     |
| `class`             | Int16   | 1    | [`nparticle`]| NA    | Orbit class of the MC particle              |
| `weight`            | Float64 | 1    | [`nparticle`]| fast-ions| Weight of the MC particle                |
| `r`                 | Float64 | 1    | [`nparticle`]| cm    | R positions of the MC particle              |
| `z`                 | Float64 | 1    | [`nparticle`]| cm    | Z positions of the MC particle              |
| `phi`               | Float64 | 1    | [`nparticle`]| rad   | Phi positions of the MC particle (Optional) |
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
