title: FIDASIM Outputs

This page documents FIDASIMs Outputs. It is under constructions

[TOC]

---
# NPA Outputs
The npa output file is: `<runID>_npa.h5`. This file is only created if `n_npa > 0` or if `n_pnpa > 0`

It contains:

|       Variable      | Present if |  Type   | Rank | Dimensions | Units |               Description                |
|:-------------------:|:----------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `count`             | `n_npa > 0`| Int32   | 1    | NA         | NA    | Number of active markers that hit the detector: count(chan)|
| `pcount`            | `n_pnpa > 0`| Int32   | 1    | NA         | NA    | Number of passive markers that hit the detector: count(chan)|
| `energy`            | `n_npa > 0`| Float64 | 1    | nE         | keV   | Energy axis where the NPA flux is binned |
| `flux`              | `n_npa > 0`| Float64 | 2    | nE, nchan  | #/(s*keV)| Active Neutral flux: flux(energy,chan)   |
| `pflux`             | `n_pnpa > 0`| Float64 | 2    | nE, nchan  | #/(s*keV)| Passive Neutral flux: flux(energy,chan)   |
| `nchan `            | `n_(p)npa > 0`| Int32   | 0    | NA         | NA    | Number of channels     |
| `nenergy `          | `n_(p)npa > 0`| Int32   | 0    | NA         | NA    | Number of energies     |
| `radius`            | `n_(p)npa > 0`| Float64 | 1    | nchan      | m     | Detector LOS radius at midplane or tangency point     |
| `particles`         | `n_npa > 1`|  Structure| NA   | NA         | NA    | Data from the active markers, see below           |
| `passive_particles` | `n_pnpa > 1`| Structure   | NA    | NA         | NA    | Data from the passive markers, see below      |

If the NPA switch is greater than 1, the `particle` or `passive_particle` structure is also output in the file, it contains, in both cases:

|       Variable      | Present if    |  Type   | Rank | Dimensions | Units |               Description                |
|:-------------------:|:-------------:|:-------:|:----:|:----------:|:-----:|:-----------------------------------------|
| `nparticle`         | `n_(p)npa > 1`| Int32   | 0    | NA         | NA    | Number of markers that hit the detector|
| `class`             | `n_(p)npa > 1`| Int32   | 1    | nMarker    | NA    | Class of the neutral particle |
| `detector`          | `n_(p)npa > 1`| Int32   | 1    | nMarker    | NA    | Detector (channel) the marker hit|
| `energy`            | `n_(p)npa > 1`| Float64 | 1    | nMarker    | keV   | Energy of the particle |
| `gci`               | `n_(p)npa > 1`| Float64 | 2    | 3, nMarker | cm     | Neutral particle's gyrocenter birth position in machine coordinates: gci([x,y,z],particle)|
| `ri`                | `n_(p)npa > 1`| Float64 | 2    | 3, nMarker | cm     | Neutral particle's birth position in machine coordinates: ri([x,y,z],particle)|
| `rf`                | `n_(p)npa > 1`| Float64 | 2    | 3, nMarker | cm     | Neutral particle's hit position in machine coordinates: rf([x,y,z],particle)|
| `pitch`             | `n_(p)npa > 1`| Float64 | 1    | nMarker    | NA    | Pitch value of the neutral particle: p = v_parallel/v w.r.t. the magnetic field   |
| `weight`            | `n_(p)npa > 1`| Float64 | 1    | nMarker    | #/s   | Neutral particle's contribution to the flux |

