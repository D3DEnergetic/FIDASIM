title: PREFIDA Inputs

#PREFIDA Inputs

[TOC]

#Inputs Structure
The `inputs` structure contains the basic information needed by FIDASIM.

##General Settings
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

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `shot`              | Int32   | 0    | NA         | NA    | Shot Number                              |
| `time`              | Float64 | 0    | NA         | s     | Time                                     |
| `runid`             | String  | 0    | NA         | NA    | Run ID                                   |
| `comment`           | String  | 0    | NA         | NA    | Comment                                  |
| `result_dir`        | String  | 0    | NA         | NA    | Result directory                         |
| `tables_file`       | String  | 0    | NA         | NA    | Atomic Tables file                       |

##Simulation Switches
The simulation switches can take values 0, 1, or 2. A value of zero and one will turn the calculation off and on respectively. 
A value of two will turn on additional functionality.

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `calc_bes`          | Int16   | 0    | NA         | NA    | Calculate Beam Emission and Halo spectra |
| `calc_brems`        | Int16   | 0    | NA         | NA    | Calculate Bremsstrahlung                 |
| `calc_fida`         | Int16   | 0    | NA         | NA    | Calculate FIDA spectra                   |
| `calc_npa`          | Int16   | 0    | NA         | NA    | Calculate NPA flux                       |
| `calc_birth`        | Int16   | 0    | NA         | NA    | Calculate Birth profile                  |
| `calc_fida_wght`    | Int16   | 0    | NA         | NA    | Calculate FIDA weight functions          | 
| `calc_npa_wght`     | Int16   | 0    | NA         | NA    | Calculate NPA weight functions           |
| `dump_dcx`          | Int16   | 0    | NA         | NA    | Dump DCX neutrals and spectra            |

##Monte Carlo Settings
These settings control the number of Monte Carlo particles used by FIDASIM.
Using too few particles will execute quickly but will be extremely noisy.
Contrarily, using too many particles will increase runtime but will have small Monte Carlo noise.
The following settings provide a good balance between runtime and Monte Carlo noise.

* `n_fida` = 5000000L
* `n_npa` = 5000000L
* `n_nbi` = 50000L
* `n_dcx` = 500000L
* `n_halo` = 500000L
* `n_birth` = 50000L

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `n_fida`            | Int32   | 0    | NA         | NA    | Number of FIDA MC particles              |
| `n_npa`             | Int32   | 0    | NA         | NA    | Number of NPA MC particles               |
| `n_nbi`             | Int32   | 0    | NA         | NA    | Number of NBI MC particles               |
| `n_dcx`             | Int32   | 0    | NA         | NA    | Number of DCX MC particles               |
| `n_halo`            | Int32   | 0    | NA         | NA    | Number of HALO MC particles              |
| `n_birth`           | Int32   | 0    | NA         | NA    | Number of Birth particles outputed       |

##Neutral Beam Settings
These variables define the neutral beam properties.
Currently the mass of the beam species, `ab`,  can only be the mass either protium or deuterium.
The `current_fractions` variable must sum to one.

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `ab`                | Float64 | 0    | NA         | amu   | Beam species mass                        |
| `pinj`              | Float64 | 0    | NA         | MW    | Beam power                               |
| `einj`              | Float64 | 0    | NA         | keV   | Beam energy                              |
| `current_fractions` | Float64 | 1    | [3]        | NA    | Current fractions (Full, Half, Third)    |

##Plasma Settings
These variables define the properties of the thermal plasma species.
Like the `ab` variable, `ai` can only be the mass of either protium or deuterium

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `ai`                | Float64 | 0    | NA         | amu   | Thermal Ion species mass                 |
| `impurity_charge`   | Int16   | 0    | NA         | NA    | Impurity Charge                          |

##Beam Grid Settings

These variables define a rotated coordinate system. 
The [Tait-Bryan rotation angles](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) (`alpha`,`beta`,`gamma`) define a intrinsic rotation matrix, \(R\), that is used to transform from beam grid coordinates(xyz) to machine coordinates(uvw)
$$ \mathrm{uvw = R \cdot xyz + origin} $$
If the rotation angles and `origin` are set to zero then the rotation matrix is the Identity matrix and the coordinate system is identical to machine coordinates. 


Understanding these variables can be difficult and can best be described by an example. 

1. With your right hand point your index finger pointing in the +x direction with your middle finger and thumb pointing in the +y and +z direction respectively.
2. Rotate about your thumb (z-axis) by `alpha` (ccw = +angle, cw = -angle)
3. Rotate about your middle finger (y'-axis) by `beta`
4. Rotate about your index finger (x"-axis) by `gamma`
5. Move your right hand to the `origin`
6. Define `(x|y|z)_(min|max)` by this coordinate system with your index finger being the new +x-axis

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `nx`                | Int16   | 0    | NA         | NA    | Number of cells in the X direction       |
| `ny`                | Int16   | 0    | NA         | NA    | Number of cells in the Y direction       |
| `nz`                | Int16   | 0    | NA         | NA    | Number of cells in the Z direction       |
| `xmin`              | Float64 | 0    | NA         | cm    | Minimum X value in beam grid coordinates |
| `xmax`              | Float64 | 0    | NA         | cm    | Maximum X value in beam grid coordinates |
| `ymin`              | Float64 | 0    | NA         | cm    | Minimum Y value in beam grid coordinates |
| `ymax`              | Float64 | 0    | NA         | cm    | Maximum Y value in beam grid coordinates |
| `zmin`              | Float64 | 0    | NA         | cm    | Minimum Z value in beam grid coordinates |
| `zmax`              | Float64 | 0    | NA         | cm    | Maximum Z value in beam grid coordinates |
| `alpha`             | Float64 | 0    | NA         | rad   | Tait-Bryan rotation angle about z-axis   |
| `beta`              | Float64 | 0    | NA         | rad   | Tait-Bryan rotation angle about y'-axis  |
| `gamma`             | Float64 | 0    | NA         | rad   | Tait-Bryan rotation angle about x"-axis  |
| `origin`            | Float64 | 1    | [3]        | cm    | Beam grid origin in Machine Coordinates  |

##Wavelength Grid Settings
These variables define the wavelength grid. Using a fine wavelength has no performance penalty.

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `nlambda`           | Int16   | 0    | NA         | NA    | Number of wavelengths                    |
| `lambdamin`         | Float64 | 0    | NA         | nm    | Minimum wavelength                       |
| `lambdamax`         | Float64 | 0    | NA         | nm    | Maximum wavelength                       |

##Weight Function Settings
|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `ne_wght`           | Int16   | 0    | NA         | NA    | Number of weight function energies       |
| `np_wght`           | Int16   | 0    | NA         | NA    | Number of weight function pitches        |
| `nphi_wght`         | Int16   | 0    | NA         | NA    | Number of gyro-angles                    |
| `emax_wght`         | Float64 | 0    | NA         | keV   | Maximum energy of weight functions       |
| `nlambda_wght`      | Int16   | 0    | NA         | NA    | Number of weight function wavelengths    |
| `lambdamin_wght`    | Float64 | 0    | NA         | nm    | Minimum weight function wavelength       |
| `lambdamax_wght`    | Float64 | 0    | NA         | nm    | Maximum weight function wavelength       |

#Interpolation Grid Structure
The `grid` structure contains the definition of the 2D R-Z grid that the plasma parameters and electromagnetic fields are mapped onto. 

|       Variable      |   Type  | Rank |  Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:-----------------------------------------|
| `nr`                | Int16   | 0    | NA          | NA    | Number of radii                          |
| `nz`                | Int16   | 0    | NA          | NA    | Number of z values                       |
| `r`                 | Float64 | 1    | [`nr`]      | cm    | Array of radii                           |
| `z`                 | Float64 | 1    | [`nz`]      | cm    | Array of z values                        |
| `r2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of radii `r = r2d(r,z)`         |
| `z2d`               | Float64 | 2    | [`nr`,`nz`] | cm    | 2D array of z values `z = z2d(r,z)`      |


#Neutral Beam Geometry Structure
The `nbi` structure contains the neutral beam geometry.
The `shape` of the source grid take the value of 1 or 2 for a rectangular and circular source grid respectively.

|       Variable      |   Type  | Rank |  Dimensions | Units |                    Description                     |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:---------------------------------------------------|
| `name`              | String  | 0    | NA          | NA    | Name of the neutral beam                           |
| `shape`             | Int16   | 0    | NA          | NA    | Shape of the beam source grid (1 or 2)             |
| `data_source`       | String  | 0    | NA          | NA    | Source of the neutral beam geometry                |
| `src`               | Float64 | 1    | [3]         | cm    | Position of the source grid in machine coordinates |
| `axis`              | Float64 | 1    | [3]         | NA    | Direction of the beam center line                  |
| `widy`              | Float64 | 0    | NA          | cm    | Source grid half-width in the horizontal direction |
| `widz`              | Float64 | 0    | NA          | cm    | Source grid half-height in the vertical direction  |
| `divy`              | Float64 | 1    | [3]         | rad   | Horizontal beam divergence                         |
| `divz`              | Float64 | 1    | [3]         | rad   | Vertical beam divergence                           |
| `focy`              | Float64 | 0    | NA          | cm    | Horizontal focal length                            |
| `focz`              | Float64 | 0    | NA          | cm    | Vertical focal length                              |


#Fields Structure
This structure contain the electromagnetic fields mapped onto the interpolation grid.

|       Variable      |   Type  | Rank |  Dimensions | Units |                           Description                         |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:--------------------------------------------------------------| 
| `time`              | Float64 | 0    | NA          | s     | Time when the fields data were collected/reconstructed        |
| `data_source`       | String  | 0    | NA          | NA    | Source of the fields data                                     |
| `mask`              | Int16   | 2    | [`nr`,`nz`] | NA    | Boolean mask that indicates where the fields are well defined |
| `br`                | Float64 | 2    | [`nr`,`nz`] | T     | Radial component of the magnetic field                        |
| `bt`                | Float64 | 2    | [`nr`,`nz`] | T     | Torodial/Phi component of the magnetic field                  |
| `bz`                | Float64 | 2    | [`nr`,`nz`] | T     | Z component of the magnetic field                             |
| `er`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Radial component of the electric field                        |
| `et`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Torodial/Phi component of the electric field                  |
| `ez`                | Float64 | 2    | [`nr`,`nz`] | V/m   | Z component of the electric field                             |

#Plasma Structure
This structure contain the plasma parameters mapped onto the interpolation grid.

|       Variable      |   Type  | Rank |  Dimensions | Units |                           Description                         |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:--------------------------------------------------------------| 
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

#Distribution Structure
The `dist` structure contains the fast-ion distribution which can be one of three different types.

##Fast-ion Distribution Function
|       Variable      |   Type  | Rank |           Dimensions           |          Units         |           Description           |
|:-------------------:|:-------:|:----:|:------------------------------:|:----------------------:|:--------------------------------| 
| `type`              | Int16   | 0    | NA                             | NA                     | Distribution type (1)           |
| `time`              | Float64 | 0    | NA                             | s                      | Time of the distribution        |
| `data_source`       | String  | 0    | NA                             | NA                     | Source of the distribution data |
| `nr`                | Int16   | 0    | NA                             | NA                     | Number of radii                 |
| `nz`                | Int16   | 0    | NA                             | NA                     | Number of z values              |
| `nenergy`           | Int16   | 0    | NA                             | NA                     | Number of energy values         |
| `npitch`            | Int16   | 0    | NA                             | NA                     | Number of pitch values          |
| `R`                 | Float64 | 1    | [`nr`]                         | cm                     | R array                         |
| `Z`                 | Float64 | 1    | [`nz`]                         | cm                     | Z array                         |
| `energy`            | Float64 | 1    | [`nenergy`]                    | keV                    | Energy array                    |
| `pitch`             | Float64 | 1    | [`npitch`]                     | NA                     | Pitch array w.r.t magnetic field|
| `denf`              | Float64 | 2    | [`nr`,`nz`]                    | cm^-3                  | Fast-ion density                |
| `f`                 | Float64 | 4    | [`nenergy`,`npitch`,`nr`,`nz`] | fast-ions/(dE dP cm^3) | Fast-ion distribution F(E,p,R,Z)|

##Guiding Center Monte Carlo Distribution
The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles.

The `class` variable can take values in the range of 1:`nclass`.


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

##Full-Orbit Monte Carlo Distribution
The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles

The `class` variable can take values in the range of 1:`nclass`.

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
| `vr`                | Float64 | 1    | [`nparticle`]| cm/s  | Radial component of the MC particle velocity |
| `vt`                | Float64 | 1    | [`nparticle`]| cm/s  | Torodial/Phi component of the MC particle velocity |
| `vz`                | Float64 | 1    | [`nparticle`]| cm/s  | Z component of the MC particle velocity |

#Spectral Geometry Structure
This structure contains the geometry of the spectroscopic systems 

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `nchan`             | Int32   | 0    | NA           | NA    | Number of channels              |
| `system`            | String  | 0    | NA           | NA    | Name of the spectrocopic system(s) |
| `data_source`       | String  | 0    | NA           | NA    | Source of the spectral geometry data |
| `id`                | String  | 1    | [`nchan`]    | NA    | Channel ID                      |
| `radius`            | Float64 | 1    | [`nchan`]    | cm    | Line of sight radius at midplane or tangency point |
| `lens`              | Float64 | 2    | [3,`nchan`]  | cm    | Lens location in machine coordinates |
| `axis`              | Float64 | 2    | [3,`nchan`]  | NA    | Optical axis/direction of the lines of sight |
| `spot_size`         | Float64 | 1    | [`nchan`]    | cm    | Radius of the collecting volume |
| `sigma_pi`          | Float64 | 1    | [`nchan`]    | NA    | Ratio of the intensities of the sigma and pi stark lines |

#NPA Geometry Structure
This structure contains the geometry of the spectroscopic systems 
The shapes of the detector and aperture can take the value 1 or 2 for a rectangular and circular aperture/detector respectively.

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `nchan`             | Int32   | 0    | NA           | NA    | Number of channels              |
| `system`            | String  | 0    | NA           | NA    | Name of the NPA system(s)       |
| `data_source`       | String  | 0    | NA           | NA    | Source of the NPA geometry data |
| `id`                | String  | 1    | [`nchan`]    | NA    | Channel ID                      |
| `radius`            | Float64 | 1    | [`nchan`]    | cm    | Line of sight radius at midplane or tangency point |
| `a_shape`           | Int16   | 1    | [`nchan`]    | NA    | Shape of the aperture           |
| `d_shape`           | Int16   | 1    | [`nchan`]    | NA    | Shape of the detector           |
| `a_cent`            | Float64 | 2    | [3,`nchan`]  | cm    | Position of the center of the aperture |
| `a_redge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the apertures right edge |
| `a_tedge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the apertures top edge |
| `d_cent`            | Float64 | 2    | [3,`nchan`]  | cm    | Position of the center of the detector |
| `d_redge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the detectors right edge |
| `d_tedge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the detectors top edge |
