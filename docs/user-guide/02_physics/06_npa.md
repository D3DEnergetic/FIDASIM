title: Neutral Particle Analyzer

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

# NPA Geometry
![Detector](|media|/npa_detector.png){: width="250"} ![Aperture](|media|/npa_aperture.png){: width="250"}
{: style="text-align: center"}

An NPA detector is defined by an aperture for which neutral particles must pass through and an detector.
The aperture/detectors are defined by three points and a shape as shown above.
It is assumed that between the aperture and the detector that particles travel in straight lines i.e. there is no stripping foil.

The full definition of the NPA detector is given below

|       Variable      |   Type  | Rank |  Dimensions  | Units |           Description           |
|:-------------------:|:-------:|:----:|:------------:|:-----:|:--------------------------------| 
| `nchan`             | Int32   | 0    | NA           | NA    | Number of channels              |
| `system`            | String  | 0    | NA           | NA    | Name of the NPA system(s)       |
| `data_source`       | String  | 0    | NA           | NA    | Source of the NPA geometry data |
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

# Monte Carlo NPA calculation 
The Monte Carlo method of calculating the NPA flux (MC-NPA) is as follows

1. Sample Fast-ion distribution function and get initial position and velocity
2. Determine if particle would go through the NPA aperture and hit the detecting region. If not increase counter and goto 1 else goto 3
3. Charge exchange the ion (set initial state) and solve the collisional radiative model along particle track.
4. Sum the final state of the neutral and bin the particle by its energy.
5. Repeat N times

This process works but is not very efficient. The main issue is that more often or not a particle would not hit the detector.
This causes very bad Monte Carlo noise in the calculated NPA flux.
The only way to get around this is to use a lot of particles which can be prohibitively expensive.

It would be much faster to just fire the particles directly at the NPA detector and then scale the resultant flux by the probability of that trajectory occuring.
This is approach taken in the weight function method (WF-NPA) detailed [here](./07_weights.html#npa).

An example of the calculated NPA flux for the two different methods are shown below.

![NPA Flux](|media|/npa.png)
{: style="text-align: center"}

As you can see the WF-NPA method produces superior results at a fraction of the runtime.

# Relevant Namelist Settings
* `n_npa`: Number of Monte Carlo particles used in MC-NPA calculation 
* `calc_npa`: Calculate NPA flux using the Monte Carlo Method
* `calc_npa_wght`: Calculate NPA weight function and flux using the weight function method
* `ne_wght`: Number of energies in weight function calculation
* `np_wght`: Number of pitches in weight function calculation
* `emax_wght`: Maximum energy in weight function calculation

# Fortran References
* [[read_npa]]: Reads NPA geometry and calculates NPA geometric factor
* [[npa_f]]: MC-NPA routine
* [[npa_weights]]: WF-NPA routine
