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
| `calc_bes`          | Int16   | 0    | NA         | NA    | Calculate NBI Spectra                    |
| `calc_dcx`          | Int16   | 0    | NA         | NA    | Calculate Direct Charge Exchange Spectra |
| `calc_halo`         | Int16   | 0    | NA         | NA    | Calculate HALO spectra                   |
| `calc_cold`         | Int16   | 0    | NA         | NA    | Calculate COLD spectra                   |
| `calc_brems`        | Int16   | 0    | NA         | NA    | Calculate Bremsstrahlung                 |
| `calc_fida`         | Int16   | 0    | NA         | NA    | Calculate FIDA spectra                   |
| `calc_npa`          | Int16   | 0    | NA         | NA    | Calculate NPA flux                       |
| `calc_pfida`        | Int16   | 0    | NA         | NA    | Calculate passive FIDA spectra           |
| `calc_pnpa`         | Int16   | 0    | NA         | NA    | Calculate passive NPA flux               |
| `calc_neutron`      | Int16   | 0    | NA         | NA    | Calculate B-T Neutron rate               |
| `calc_birth`        | Int16   | 0    | NA         | NA    | Calculate Birth profile                  |
| `calc_fida_wght`    | Int16   | 0    | NA         | NA    | Calculate FIDA weight functions          |
| `calc_npa_wght`     | Int16   | 0    | NA         | NA    | Calculate NPA weight functions           |

##Monte Carlo Settings
These settings control the number of Monte Carlo particles used by FIDASIM.
Using too few particles will execute quickly but will be extremely noisy.
Contrarily, using too many particles will increase runtime but will have small Monte Carlo noise.
The following settings provide a good balance between runtime and Monte Carlo noise.

* `n_fida` = 5000000L
* `n_pfida` = 50000000L
* `n_npa` = 5000000L
* `n_pnpa` = 50000000L
* `n_nbi` = 50000L
* `n_halo` = 500000L
* `n_dcx` = 500000L
* `n_birth` = 10000L

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `n_fida`            | Int64   | 0    | NA         | NA    | Number of FIDA MC particles              |
| `n_pfida`           | Int64   | 0    | NA         | NA    | Number of passive FIDA MC particles      |
| `n_npa`             | Int64   | 0    | NA         | NA    | Number of NPA MC particles               |
| `n_pnpa`            | Int64   | 0    | NA         | NA    | Number of passive NPA MC particles       |
| `n_nbi`             | Int64   | 0    | NA         | NA    | Number of NBI MC particles               |
| `n_halo`            | Int64   | 0    | NA         | NA    | Number of HALO MC particles              |
| `n_dcx`             | Int64   | 0    | NA         | NA    | Number of DCX MC particles               |
| `n_birth`           | Int64   | 0    | NA         | NA    | Number of Birth particles outputed       |
###NPA switch
* `n_npa == 0`: Will deactivate the NPA calculation
* `n_npa == 1`: Will perform the NPA calculation but just store as output the flux (per unit of energy) of neutral particle arriving to the detector
* `n_npa == 2`: Will do the same of `n_npa == 1` but will also save the initial position and velocity (Energy and pitch) of each neutral marker. Since Nov 2024, the position and energy of the gyrocenter is also saved, if the used fast-ion distribution as input is symmetric in toroidal angle
##Neutral Beam Settings
These variables define the neutral beam properties.
Currently the mass of the beam species, `ab`, can only be the mass either protium or deuterium.
The `current_fractions` variable must sum to one.
Click [here](../02_physics/04_neutrals.html#neutral-beam-density) for more information.

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `ab`                | Float64 | 0    | NA         | amu   | Beam species mass                        |
| `pinj`              | Float64 | 0    | NA         | MW    | Beam power                               |
| `einj`              | Float64 | 0    | NA         | keV   | Beam energy                              |
| `current_fractions` | Float64 | 1    | [3]        | NA    | Current fractions (Full, Half, Third)    |

##Beam Grid Settings

These variables define a rotated coordinate system.
Click [here](../02_physics/04_neutrals.html#beam-grid) for more details.

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
| `beta`              | Float64 | 0    | NA         | rad   | Tait-Bryan rotation angle about yp-axis  |
| `gamma`             | Float64 | 0    | NA         | rad   | Tait-Bryan rotation angle about xpp-axis |
| `origin`            | Float64 | 1    | [3]        | cm    | Beam grid origin in Machine Coordinates  |

##Wavelength Grid Settings
These variables define the wavelength grid. Using a fine wavelength has no performance penalty.
Click [here](../02_physics/05_spectra.html) for more more information.

|       Variable      |   Type  | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `nlambda`           | Int16   | 0    | NA         | NA    | Number of wavelengths                    |
| `lambdamin`         | Float64 | 0    | NA         | nm    | Minimum wavelength                       |
| `lambdamax`         | Float64 | 0    | NA         | nm    | Maximum wavelength                       |

##Weight Function Settings
These variables define the setting for the calculation of weight functions.
Click [here](../02_physics/07_weights.html) for more information.

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
The `grid` structure contains the definition of the 2D/3D cylindrical grid that the plasma parameters and electromagnetic fields are mapped onto.

**Note:** If `nphi` is specified without the `phi` array, the grid remains 2D for plasma/field interpolation, but `nphi` is used to set the toroidal resolution of the passive grid. This allows finer toroidal resolution for passive diagnostics without the computational overhead of full 3D plasma calculations.

|       Variable      |   Type  | Rank |  Dimensions | Units |               Description                |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:-----------------------------------------|
| `nr`                | Int16   | 0    | NA          | NA    | Number of radii                          |
| `nz`                | Int16   | 0    | NA          | NA    | Number of z values                       |
| `nphi`              | Int16   | 0    | NA          | NA    | Number of phi values (Optional)          |
| `r`                 | Float64 | 1    | [`nr`]      | cm    | Array of radii                           |
| `z`                 | Float64 | 1    | [`nz`]      | cm    | Array of z values                        |
| `phi`               | Float64 | 1    | [`nphi`]    | rad   | Array of phi values (Optional)           |


#Neutral Beam Geometry Structure
The `nbi` structure contains the neutral beam geometry.
The `(a)shape` of the source grid and apertures take the value of 1 or 2 for a rectangular and circular respectively.
Click [here](../02_physics/04_neutrals.html#neutral-beam-geometry) for more information.

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
| `naperture`         | Int16   | 0    | NA          | NA    | Number of apertures                                |
| `ashape`            | Int16   | 1    |[`naperture`]| NA    | Shape of the aperture(s) (1 or 2)                  |
| `awidy`             | Float64 | 1    |[`naperture`]| cm    | Half-width of the aperture(s)                      |
| `awidz`             | Float64 | 1    |[`naperture`]| cm    | Half-height of the aperture(s)                     |
| `aoffy` | Float64 | 1 |[`naperture`]| cm | Horizontal (y) offset of the aperture(s) relative to the +x aligned beam centerline |
| `aoffz` | Float64 | 1 |[`naperture`]| cm | Vertical (z) offset of the aperture(s) relative to the +x aligned beam centerline |
| `adist` | Float64 | 1 |[`naperture`]| cm | Distance from the center of the beam source grid to the aperture(s) plane |

#Fields Structure
This structure contain the electromagnetic fields mapped onto the interpolation grid.
Click [here](../02_physics/01_plasma.html#plasma-parameters-and-fields) for more information.

|       Variable      |   Type  | Rank |  Dimensions | Units |                           Description                         |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:--------------------------------------------------------------| 
| `time`              | Float64 | 0    | NA          | s     | Time when the fields data were collected/reconstructed        |
| `data_source`       | String  | 0    | NA          | NA    | Source of the fields data                                     |
| `mask`              | Int16   | 2/3    | [`nr`,`nz`[,`nphi`]] | NA    | Boolean mask that indicates where the fields are well defined |
| `br`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | T     | Radial component of the magnetic field                        |
| `bt`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | T     | Torodial/Phi component of the magnetic field                  |
| `bz`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | T     | Z component of the magnetic field                             |
| `er`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | V/m   | Radial component of the electric field                        |
| `et`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | V/m   | Torodial/Phi component of the electric field                  |
| `ez`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | V/m   | Z component of the electric field                             |
| `description`       | String  | 0    | NA          | NA    | Electromagnetic Field                                         |
| `coordinate system` | String  | 0    | NA          | NA    | Cylindrical                                                   |

#Plasma Structure
This structure contain the plasma parameters mapped onto the interpolation grid.
Click [here](../02_physics/01_plasma.html#plasma-parameters-and-fields) for more information.

|       Variable      |   Type  | Rank |  Dimensions | Units |                           Description                         |
|:-------------------:|:-------:|:----:|:-----------:|:-----:|:--------------------------------------------------------------| 
| `time`              | Float64 | 0    | NA          | s     | Time when the plasma parameter data was collected             |
| `data_source`       | String  | 0    | NA          | NA    | Source of the plasma parameter data                           |
| `mask`              | Int16   | 2/3    | [`nr`,`nz`[,`nphi`]] | NA    | Boolean mask that indicates where the plasma is well defined  |
| `te`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | keV   | Electron temperature                                          |
| `ti`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | keV   | Ion temperature                                               |
| `dene`              | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | cm^-3 | Electron density                                              |
| `zeff`              | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | NA    | Z-effective                                                   |
| `vr`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | cm/s  | Radial component of the bulk plasma rotation/flow             |
| `vt`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | cm/s  | Torodial/Phi component of the bulk plasma rotation/flow       |
| `vz`                | Float64 | 2/3    | [`nr`,`nz`[,`nphi`]] | cm/s  | Z component of the bulk plasma rotation/flow                  |
| `description`       | String  | 0    | NA          | NA    | Plasma Parameters                                             |
| `coordinate system` | String  | 0    | NA          | NA    | Cylindrical                                                   |

#Distribution Structure
The `dist` structure contains the fast-ion distribution which can be one of three different types.
Click [here](../02_physics/01_plasma.html#distributions) for more information.

##Fast-ion Distribution Function
|       Variable      |   Type  | Rank |           Dimensions           |          Units         |           Description           |
|:-------------------:|:-------:|:----:|:------------------------------:|:----------------------:|:--------------------------------| 
| `type`              | Int16   | 0    | NA                                      | NA                     | Distribution type (1)                  |
| `r`                 | Float64 | 1    | [`nr`]                                  | cm                     | Array of radii                         |
| `z`                 | Float64 | 1    | [`nz`]                                  | cm                     | Array of z values                      |
| `phi`               | Float64 | 1    | [`nphi`]                                | cm                     | Array of phi values (Optional)         |
| `time`              | Float64 | 0    | NA                                      | s                      | Time of the distribution               |
| `data_source`       | String  | 0    | NA                                      | NA                     | Source of the distribution data        |
| `nenergy`           | Int16   | 0    | NA                                      | NA                     | Number of energy values                |
| `npitch`            | Int16   | 0    | NA                                      | NA                     | Number of pitch values                 |
| `energy`            | Float64 | 1    | [`nenergy`]                             | keV                    | Energy array                           |
| `pitch`             | Float64 | 1    | [`npitch`]                              | NA                     | Pitch array w.r.t magnetic field       |
| `denf`              | Float64 | 3    | [`nr`,`nz`[,`nphi`]]                    | cm^-3                  | Fast-ion density                       |
| `f`                 | Float64 | 5    | [`nenergy`,`npitch`,`nr`,`nz`[,`nphi`]] | fast-ions/(dE dP cm^3) | Fast-ion distribution F(E,p,R,Z[,Phi]) |

##Guiding Center Monte Carlo Distribution
The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles.

The `class` variable can take values in the range of 1:`nclass`.


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

##Full-Orbit Monte Carlo Distribution
The sum(`weight`) = # of Fast-ions in phase space sampled by the MC particles

The `class` variable can take values in the range of 1:`nclass`.

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `type`              | Int16   | 0    | NA           | NA    | Distribution type (3)                              |
| `time`              | Float64 | 0    | NA           | s     | Time of the distribution                           |
| `data_source`       | String  | 0    | NA           | NA    | Source of the distribution data                    |
| `nparticle`         | Int32   | 0    | NA           | NA    | Number of MC particles                             |
| `nclass`            | Int16   | 0    | NA           | NA    | Number of orbit classes                            |
| `class`             | Int16   | 1    | [`nparticle`]| NA    | Orbit class of the MC particle                     |
| `weight`            | Float64 | 1    | [`nparticle`]| fast-ions| Weight of the MC particle                       |
| `r`                 | Float64 | 1    | [`nparticle`]| cm    | R positions of the MC particle                     |
| `z`                 | Float64 | 1    | [`nparticle`]| cm    | Z positions of the MC particle                     |
| `phi`               | Float64 | 1    | [`nparticle`]| rad   | Phi positions of the MC particle (Optional)        |
| `vr`                | Float64 | 1    | [`nparticle`]| cm/s  | Radial component of the MC particle velocity       |
| `vt`                | Float64 | 1    | [`nparticle`]| cm/s  | Torodial/Phi component of the MC particle velocity |
| `vz`                | Float64 | 1    | [`nparticle`]| cm/s  | Z component of the MC particle velocity            |

#Spectral Geometry Structure
This structure contains the geometry of the spectroscopic systems 
Click [here](../02_physics/05_spectra.html) for more more information.

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
Click [here](../02_physics/06_npa.html) for more more information.

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

#Adaptive Time Step Settings
These variables define the calculations of an adaptive time step in the track and track_cylindrical subroutines.
The number of splits in a cell is determined by the equation:
$$ n_{cells} = ceil(\frac{\Delta x}{x_{avg} \cdot dl \cdot tol}) $$
The value of \(dl\) is calculated by \(dt \cdot vn\), dt is the time step and vn is the normal component of velocity.
The tol corresponds to `split_tol`, defined as the percent change per cm.
The variable x is a plasma parameter defined by `adaptive` according to the following values:

* 0: Adaptive off
* 1: `dene`, electron density
* 2: `denn`, cold neutral density averaged over ground states of thermal species
* 3: `denf`, fast-ion density
* 4: `deni`, ion density averaged over thermal species
* 5: `denimp`, impurity density
* 6: `te`, electron temperature
* 7: `ti`, ion temperature

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
|`adaptive`           |Int32    |0     |NA            |NA     |Calculate `n_cells` according to plasma parameter|
|`split_tol`          |Float64  |0     |NA            |cm^-1  |Split tolerance, fractional change/cm|
|`max_cell_splits`    |Int32    |0     |NA            |NA     |Upper limit for n_cells          |
