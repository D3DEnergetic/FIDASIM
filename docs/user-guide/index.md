title: User Guide

#Introduction
In fusion plasmas fast-ions can undergo the following process with injected neutral Hydrogen:
$$ H^+ + H(n) \rightarrow H^\ast(m) + H^+ $$
where \(H^\ast(m)\) is an excited state of Hydrogen.
The newly created fast-neutral \(H^\ast(m)\) can then be collisionally de-excited to a lower state, emitting a doppler shifted photon as illustrated below.

![Charge Exchange Process](|media|/fast_ion_process.svg Charge Exchange Process){: width="500" }
{: style="text-align: center" }

Both the fast-neutral and the photon contain information about the velocity of fast-ion before it was neutralized.
There are two types diagnostics that take advantage of this fact: Fast-ion D-α (FIDA) diagnostics and Neutral Particle Analyzers (NPA).
The interpretation of the diagnostic signals requires modeling of the above process which has a complicated dependence on the plasma parameters, electromagnetic fields, and neutral beam and diagnostic geometry.

Towards this end, FIDASIM was created.  

##History
The very first implementation of FIDASIM was written by Yadong Luo as a part in his [thesis](http://www.physics.uci.edu/~wwheidbr/papers/thesis_luo.pdf).
The code would be later be improved upon by [Bill Heidbrink and Deyong Lui](http://www.physics.uci.edu/~wwheidbr/papers/FIDASIM.pdf) for public use.

Originally, FIDASIM was written in the IDL programming language and was prohibitively slow.
As a part of his [thesis](http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/46/051/46051941.pdf) Ben Geiger wrote a version of FIDASIM written in Fortran 90.
This prototype version was parallelized using OpenMP and was orders of magnitude faster but was not as easy to use as the IDL version and was difficult to port to different devices. 

Most recently, [Luke Stagner](http://github.com/lstagner) as a part of his thesis has rewritten Ben Geigers Fortran 90 version to be compatible with any axisymmetric fusion device as well as additional functionality.
Luke has also put special effort on making FIDASIM user friendly; the result of which you are currently reading and should be eternally grateful. 

##Capabilities
Currently, FIDASIM has routines for calculating:

* Neutral beam deposition and birth profile
* Halo neutral density
* Visible Bremsstrahlung
* Beam and Halo spectra
* FIDA spectra
* NPA Flux
* FIDA and NPA phase-space sensitivities i.e. weight functions

##Installation
For installation instructions check out our [Getting Started](./01_getting_started.html) guide.
