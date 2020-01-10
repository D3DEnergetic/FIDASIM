title: Cold Neutrals and Passive Signals

Originally, FIDASIM assumed axisymmetry and included active signals produced only by charge exchange (CX) with injected neutrals.
However, passive signals produced by CX with cold neutrals can also be important.
Therefore, FIDASIM is improved to predict passive signals for a given edge-cold neutral population.

[TOC]

---

# Introduction

Active diagnostic signals are produced from charge exchange with injected neutrals, and passive signals are produced from charge exchange with cold neutrals.
Passive signals must be treated to get valid active FIDA data.
For example, passive-FIDA (p-FIDA) signals of comparable magnitude with active signals were experimentally measured on NSTX-U.

![NSTX-U](|media|/nstxu.png){: width="400"}
{: style="text-align: center"}

# Time evolution of cold neutrals

FIDASIM accepts 2D and 3D cold neutral density input (TRANSP output variable is `dn0wd`).
FIDASIM will assume that all neutrals are in the ground state.
Using local plasma parameters, the neutrals are time evolved by solving the collisional radiative model until equilibrium is achieved.
Then, the neutrals are distributed throughout the interpolation grid for subsequent passive calculations.

# Types of Passive Signals

FIDASIM can calculate passive signals for the following diagnostics

* Fast-ion D-\(\alpha\) (FIDA)
* Neutral Particle Analyzer (NPA)

# Passive Grid

By default, FIDASIM will define the passive grid to encompass the beam grid and the entire viewable plasma volume.
If the interpolation grid is 3D, then the passive grid is the interpolation grid.
FIDASIM writes the passive grid settings are printed to the [standard output](../01_getting_started/03_running.html#running-interactively)

# Relevant Namelist Settings

* `calc_pfida`: Calculate pFIDA spectra
* `n_pfida`: Number of Monte Carlo particles used in passive FIDA spectra calculation
* `nlambda`: Number of wavelength bins
* `lambdamin`: Minimum wavelength [nm]
* `lambdamax`: Maximum wavelength [nm]

# Fortran References

* [[make_passive_grid]]: Creates passive grid from input geometries
* [[read_chords]]: Reads FIDA geometry
* [[pfida_f]]: Calculates FIDA light using a distribution function
* [[pfida_mc]]: Calculates FIDA light using a Monte Carlo fast-ion distribution
* [[read_npa]]: Reads NPA geometry and calculates NPA geometric factor
* [[pnpa_f]]: Calculates passive MC-NPA routine using a distribution function
* [[pnpa_mc]]: Calculates passive NPA flux using a Monte Carlo fast-ion distribution
