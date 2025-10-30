title: Neutron Collimator

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

[TOC]

---

# Neutron Collimator Geometry

The neutron collimator geometry is defined similarly to the NPA geometry, with an aperture and a detector.

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------|
| `nchan`             | Int32   | 0    | NA           | NA    | Number of channels              |
| `system`            | String  | 0    | NA           | NA    | Name of the NC system(s)       |
| `data_source`       | String  | 0    | NA           | NA    | Source of the NC geometry data |
| `id`                | String  | 1    | [`nchan`]    | NA    | Channel ID                      |
| `radius`            | Float64 | 1    | [`nchan`]    | cm    | Line of sight radius at midplane or tangency point |
| `a_shape`           | Int16   | 1    | [`nchan`]    | NA    | Shape of the aperture (1=rect, 2=circ) |
| `d_shape`           | Int16   | 1    | [`nchan`]    | NA    | Shape of the detector (1=rect, 2=circ) |
| `a_cent`            | Float64 | 2    | [3,`nchan`]  | cm    | Position of the center of the aperture |
| `a_redge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the apertures right edge |
| `a_tedge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the apertures top edge |
| `d_cent`            | Float64 | 2    | [3,`nchan`]  | cm    | Position of the center of the detector |
| `d_redge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the detectors right edge |
| `d_tedge`           | Float64 | 2    | [3,`nchan`]  | cm    | Position of the detectors top edge |

# Calculation

The neutron collimator calculation in FIDASIM provides two main outputs: the neutron flux and the weight functions. The calculation is enabled by setting the `calc_nc_wght` parameter in the `fidasim_inputs` namelist.

## Neutron Flux

The neutron flux is the number of neutrons per second per unit energy that reach the detector. It is calculated by integrating the neutron source over the plasma volume, weighted by the solid angle of the detector and the transmission probability of the collimator.

## Weight Functions

The neutron collimator weight function represents the sensitivity of the detector to fast ions at different energies and pitches. It is defined as the neutron flux per fast ion, and is given in units of [neutrons/(s*fast-ion*dE*dP)]. The weight function is a powerful tool for analyzing the fast-ion distribution function.

## Emissivity

If `calc_nc_wght` is set to 2, FIDASIM will also calculate the channel-resolved neutron emissivity. This is the neutron source strength per unit volume, resolved for each detector channel, and has units of [neutrons/(s*cm^3)].

# Output File

The results of the neutron collimation calculation are saved to an HDF5 file named `<runid>_nc_weights.h5` in the specified result directory. This file contains the following datasets under the `/ncweight` group:

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------|
| `energy`            | Float64 | 1    | [`ne_nc`]    | keV   | Energy grid for the weight function and flux |
| `pitch`             | Float64 | 1    | [`np_nc`]    | NA    | Pitch grid for the weight function |
| `weight`            | Float64 | 3    | [`ne_nc`,`np_nc`,`nchan`] | neutrons/(s*fast-ion*dE*dP) | NC weight function |
| `flux`              | Float64 | 2    | [`ne_nc`,`nchan`] | neutrons/(s*dE) | Neutron flux |
| `emissivity`        | Float64 | 3    | [`nr`,`nz`,`nchan`] | neutrons/(s*cm^3) | Channel-resolved neutron emissivity (if `calc_nc_wght >= 2`) |
| `r`                 | Float64 | 1    | [`nr`]       | cm    | Major radius grid for emissivity |
| `z`                 | Float64 | 1    | [`nz`]       | cm    | Vertical position grid for emissivity |

# Relevant Namelist Settings
* `calc_nc_wght`: Enable NC calculation (1 for weights and flux, 2 for emissivity as well)
* `ne_nc`: Number of energy bins for the NC weight function.
* `np_nc`: Number of pitch bins for the NC weight function.
* `emax_nc_wght`: Maximum energy for the NC weight function [keV].

# Fortran References
* [[read_neutron_collimator]]: Reads the neutron collimator geometry.
* [[neutron_weights]]: Calculates the neutron collimator weight functions.
* [[neutron_spec_f]]: Calculates the neutron collimator flux from a distribution function.
* [[neutron_spec_mc]]: Calculates the neutron collimator flux from Monte Carlo particles.
