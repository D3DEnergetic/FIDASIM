title: Preprocessing Inputs

#PREFIDA: A FIDASIM preprocessor
FIDASIM requires inputs to be in a [specified format](../../03_technical/01_inputs.html).
[PREFIDA](|url|/sourcefile/prefida.pro.html) is an IDL routine that takes the required inputs, checks them validity, and transforms them into a form FIDASIM understands.

PREFIDA is called as follows
```
IDL> prefida, inputs, grid, nbi, plasma, fields, dist, spec=spec, npa=npa
```
where arguments are defined as follows. Click the argument's description for more information.

* `inputs`: [General Settings](./01_inputs.html)
* `grid`: [Interpolation grid](./02_grid.html)
* `nbi`: [Neutral Beam Geometry](./03_nbigeo.html)
* `fields`: [Electromagnetic Fields](./04_fields.html)
* `plasma`: [Plasma Parameters](./05_plasma.html)
* `dist`: [Fast-ion Distribution](./06_distribution.html)
* `spec`: [Spectral Geometry](./07_specgeo.html)
* `npa`: [NPA Geometry](./08_npageo.html)

PREFIDA will create the following files

* [Namelist File](../../03_technical/01_inputs.html#namelist-file)
* [Geometry File](../../03_technical/01_inputs.html#geometry-file)
* [Equilibrium File](../../03_technical/01_inputs.html#equilibrium-file)
* [Distribution File](../../03_technical/01_inputs.html#distribution-file)

Most devices may have already setup helper routines to make running FIDASIM and Prefida easy. 
Click [here](../05_devices.html) to find out if someone has done your work for you.
