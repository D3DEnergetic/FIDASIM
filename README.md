[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1341369.svg)](https://doi.org/10.5281/zenodo.1341369)
[![Build Status](https://travis-ci.org/D3DEnergetic/FIDASIM.svg?branch=master)](https://travis-ci.org/D3DEnergetic/FIDASIM)

![FIDASIM](docs/media/fidasim-logo.png)

# Citing FIDASIM Source Code
Along with the FIDASIM paper, please cite the source code.
```
@article{FIDASIMpaper,
  author       = {Benedikt Geiger and Luke Stagner and William W Heidbrink and Ralph Dux and Rainer Fischer and Yutaka Fujiwara and Alvin Garcia and Asger Schou Jacobsen and Anton Jansen vanVuuren and Alexander N Karpushov and Deyong Liu and Philip Adrian Schneider and Igor Sfiligoi and Peter Zsolt Poloskei and Markus Weiland},
  title        = {Progress in modelling fast-ion D-alpha spectra and neutral particle analyzer fluxes using FIDASIM},
  journal      = {Plasma Physics and Controlled Fusion},
  url          = {http://iopscience.iop.org/10.1088/1361-6587/aba8d7},
  year         = {2020}
}
```

```
@misc{FIDASIMcode,
  author       = {Stagner, L. and Geiger, B. and Heidbrink, W.W.},
  title        = {{FIDASIM: A Neutral Beam and Fast-ion Diagnostic Modeling Suite}},
  doi          = {10.5281/zenodo.1341369},
  url          = {https://doi.org/10.5281/zenodo.1341369}
}
```

# User Documention
Click [here](http://d3denergetic.github.io/FIDASIM/) for full user documentation.

# Developer Documentation

## Contributing Guidelines
If you would like to add code to FIDASIM please follow these guidelines

1. Use spaces to indent code not tabs
2. Don't push directly to `master`. All new work should be on a feature branch prefixed with your initials, e.g. `ls/newfeature`. When you think you are done open a pull request so it can undergo code review.
3. Add inline documentation as you go. Follow existing documentation format.
4. Be as general as possible. This makes it easier to add new features down the line.

## Common Makefile Options
FIDASIM has a few different makefile options
```
make DEBUG=y #Turns off all optimizations and OpenMP. Turns on bounds checking and floating point exceptions. Default: n

make USE_OPENMP=n #Turns off OpenMP but leaves other optimizations. Default: y

make PROFILE=y #Turns on gprof profiling. Turns off OpenMP. Default: n

make USE_MPI=y #Turns on MPI parallelization. Default: n

make ARCH=<TARGET> #Compile code optimized for <TARGET> e.g. ARCH=haswell/CORE-AVX2/native. Options differ depending on compiler. Generated code may not run on different architectures. Default: n
```
Run `make help` for all the available options

## Profiling
You can profile FIDASIM by building with the `PROFILE=y` build option.
With profiling enabled FIDASIM will output a `gmon.out` file.
To turn this into a profile report run
```
gprof ./fidasim gmon.out > fidasim.profile
```
This will create a report that looks something like this
```
Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 34.38    240.02   240.02 106095199     0.00     0.00  __eigensystem_MOD_hqr2
 13.75    336.04    96.01 106095199     0.00     0.00  __eigensystem_MOD_matinv
  9.62    403.18    67.14 106095199     0.00     0.00  __libfida_MOD_get_rate_matrix
  6.39    447.82    44.64 106095199     0.00     0.00  __eigensystem_MOD_hqrvec
  4.40    478.53    30.72  5049752     0.00     0.00  __libfida_MOD_interpol2d_2d_arr
  3.66    504.10    25.56  5049752     0.00     0.00  __parallel_rng_MOD_randind_w_2
  3.48    528.41    24.31 106095199     0.00     0.00  __eigensystem_MOD_balance
  3.11    550.11    21.71 116519743     0.00     0.00  __libfida_MOD_store_fida_photons
  3.07    571.58    21.47 482016993     0.00     0.00  __libfida_MOD_in_plasma
  2.99    592.48    20.90 106095199     0.00     0.00  __eigensystem_MOD_elmhes
  1.87    605.53    13.05 636571194     0.00     0.00  __eigensystem_MOD_outerprod
  1.53    616.19    10.66  5049520     0.00     0.00  __libfida_MOD_track.constprop.9
  1.47    626.45    10.26        1    10.26   696.97  __libfida_MOD_fida_f
...
```
This can be difficult to parse so it sometimes helpful to use a profile visualization tool like [gprof2dot](https://github.com/jrfonseca/gprof2dot) or the online tool [Profiler Visualizer 2000](http://gprof.jlf-hacks.com/)

## Documentation Website
The documentation is built using [FORD](https://github.com/cmacmackin/ford).
FORD takes inline documentation, prefixed by `!+`, and the markdown files located in `docs/user_guide` and wraps them all up in a pretty website.
The documentation website is automatically updated whenever a commit is added to the master branch. This is done through Travis-CI.

To build the documentation locally run
```
make docs
```
This will build the documentation website in the `docs/html` directory.

The above command will also check the website for dead links which requires [linkchecker](https://wummel.github.io/linkchecker/) to be installed.
This functionality can be disabled by passing the `CHECK_LINKS=n` build option to the `make docs` command. 
